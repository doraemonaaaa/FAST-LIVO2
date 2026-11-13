#include <iostream>
#include "initialization/imu_time_offset/imu_time_offset.hpp"
#include "initialization/imu_time_offset/scan_matching.hpp"
#include <fstream>
#include <iomanip>
#include <optional>
#include <chrono>

namespace initialization {
namespace {
using Policy = livo::ImuTimeOffsetEstimator;
using namespace imu_time_offset;

RotationCloud SampleQuery(const RotationCloud& input) {
    RotationCloud out; out.time = input.time;
    const size_t stride = std::max(size_t(1),(input.cloud->size()+3999)/4000);
    // Original point positions/timestamps chosen ONCE, before trying offsets.
    for (size_t i = 0; i < input.cloud->size(); i += stride) {
        out.cloud->push_back((*input.cloud)[i]); out.point_times.push_back(input.point_times[i]);
    }
    return out;
}
RotationCloud ReferenceMap(const RotationCloud& before, const RotationCloud& after,
                           const RotationFit& backward, const RotationFit& forward) {
    RotationCloud map;
    auto append = [&](const RotationCloud& cloud,const M3D& R,const V3D& t) {
        for (const auto& point : *cloud.cloud) {
            const V3D p = R*point.getVector3fMap().cast<double>()+t;
            PointType q{};q.x=p.x();q.y=p.y();q.z=p.z();map.cloud->push_back(q);
        }
    };
    // fits[i] maps scan i to i-1; maps are expressed at the query's IMU frame.
    append(before,backward.rotation.transpose(),-backward.rotation.transpose()*backward.translation);
    append(after,forward.rotation,forward.translation);
    return map;
}
Json::Value ScoreJson(const livo::OffsetMatchScore& s) {
    Json::Value out;
    out["scans"] = Json::UInt64(s.scans); out["registered_scans"] = Json::UInt64(s.valid_scans);
    out["query_points"] = Json::UInt64(s.points); out["matched_points"] = Json::UInt64(s.matches);
    out["match_ratio"] = s.points ? Json::Value(s.ratio()) : Json::Value();
    out["inlier_rms_m"] = std::isfinite(s.rms()) ? Json::Value(s.rms()) : Json::Value();
    out["capped_rms_m"] = std::isfinite(s.capped_rms()) ? Json::Value(s.capped_rms()) : Json::Value();
    const double balanced = Policy::BalancedScore(s);
    out["balanced_score"] = std::isfinite(balanced) ? Json::Value(balanced) : Json::Value();
    return out;
}
}

livo::ImuTimeOffsetResult EstimateImuTimeOffset(
    const std::vector<LivoScan>& scans,
    const std::vector<ImuSample::ConstPtr>& imu,
    const OffsetConfig& param,
    const Eigen::Matrix4d& imu_lidar,
    const std::filesystem::path& diagnostics,
    double window_begin, double window_end,
    const std::vector<OffsetInitializationState>& states) {
    livo::ImuTimeOffsetResult early; early.requested = true;
    if (!std::isfinite(window_begin) || !std::isfinite(window_end) ||
        window_begin<0 || window_end<=window_begin || states.empty()) {
        early.reason="invalid_stage_two_window"; return early;
    }
    for (size_t i=0;i<states.size();++i) {
        const auto& state=states[i];
        if (!std::isfinite(state.time) || state.time>window_end ||
            (i && state.time<=states[i-1].time) || !state.rotation.allFinite() ||
            !state.velocity.allFinite() || !state.gravity.allFinite() ||
            !state.gyro_bias.allFinite() || !state.acceleration_bias.allFinite()) {
            early.reason="invalid_stage_two_state"; return early;
        }
    }
    const auto candidate_shifts = Policy::CandidateShifts(param.imu_time_offset_s);
    double maximum_shift = 0;
    for (double shift : candidate_shifts) maximum_shift = std::max(maximum_shift,std::abs(shift));
    for (size_t i = 0; i < imu.size(); ++i)
        if (!imu[i] || !std::isfinite(imu[i]->timestamp) ||
            !imu[i]->angular_velocity.allFinite() || !imu[i]->acceleration.allFinite() ||
            (i && imu[i]->timestamp <= imu[i-1]->timestamp)) {
            early.reason = "invalid_imu"; return early;
        }
    std::cout << "IMU time offset initialization begin: frozen temporary LIVO states, "
                  << Policy::kWindowSeconds << " second window; offset parameter range ["
                  << Policy::kMinimumOffset*1000 << ", " << Policy::kMaximumOffset*1000
                  << "] ms, step " << Policy::kGridStep*1000 << " ms";
    std::vector<RotationCloud> full;
    double begin = -1, last = -1;
    for (const auto& scan : scans) {
        if (window_begin>=0 && scan.begin<window_begin) continue;
        if (window_end>=0 && scan.end>window_end) break;
        if (begin < 0) begin = scan.begin;
        if (scan.end-begin > Policy::kWindowSeconds) break;
        full.push_back(FullRotationCloud(RotationClouds(scan,imu_lidar))); last = scan.end;
    }
    early.scans = full.size(); early.window_s = std::max(0.,last-begin);
    early.window_begin_s = begin; early.window_end_s = last;
    if (full.size() < 5 || imu.size() < 2) {
        early.reason = "insufficient_reference_scans"; return early;
    }
    // Build reference geometry only once with LiDAR-only motion. It is never
    // rebuilt with a candidate's IMU deskew, registered query or LIVO state.
    auto corrected = full;
    std::vector<RotationFit> fits(full.size());
    std::vector<LidarTwist> twists(full.size());
    std::vector<bool> pair_enabled(full.size(),false);
    for (size_t i = 1; i < full.size(); ++i) {
        const double dt = full[i].time-full[i-1].time;
        pair_enabled[i] = !(dt <= .025 || dt >= .2);
    }
    const auto reference_start=std::chrono::steady_clock::now();
    for (int pass = 0; pass < 3; ++pass) {
        if(param.progress)param.progress();
        auto next = RegisterRotationPairs(corrected,fits,param.workers,pair_enabled);
        size_t accepted = 0;
        for (size_t i = 1; i < full.size(); ++i) accepted += next[i].valid;
        fits = std::move(next);
        // Fits are frozen. Each task updates only its own twist and cloud;
        // all clouds finish before the next registration pass starts.
        ForEachInitializationTask(full.size(),param.workers,[&](size_t i) {
            if (i && i+1 < full.size()) {
                // Replace stale valid twists even when registration fails.
                twists[i] = EstimateTwist(fits[i],full[i].time-full[i-1].time,
                    fits[i+1],full[i+1].time-full[i].time);
            }
            corrected[i] = twists[i].valid ? DeskewRotationCloud(full[i],twists[i]) : full[i];
        });
        std::clog << "[imu time offset init] reference pass=" << pass
                       << " scans=" << full.size() << " accepted=" << accepted;
    }
    // Pair final reference geometry with transforms fitted to that geometry.
    // This refit uses each pair's own frozen guess, without chain fallback.
    auto refitted = fits;
    ForEachInitializationTask(full.size()-1,param.workers,[&](size_t job) {
        const size_t i = job+1;
        const double dt = full[i].time-full[i-1].time;
        refitted[i] = dt > .025 && dt < .2 ?
            RegisterRotation(corrected[i-1],corrected[i],fits[i]) : RotationFit{};
    });
    fits = std::move(refitted);
    const auto maps_start=std::chrono::steady_clock::now();
    if(param.progress)param.progress();
    struct Query {
        size_t index; RotationCloud raw; FrozenPlaneMap map; DeskewState state; bool validation;
    };
    std::vector<Query> queries;
    const double split = begin+.7*(last-begin);
    // Even scans are queries; only odd scans enter reference maps. Use every
    // eligible query so the held-out tail retains enough observations.
    std::vector<std::optional<Query>> prepared_queries((full.size()-3)/2);
    ForEachInitializationTask(prepared_queries.size(),param.workers,[&](size_t job) {
        const size_t i = 2+2*job;
        if (!fits[i].valid || !fits[i+1].valid || !twists[i].valid ||
            !twists[i-1].valid || !twists[i+1].valid) return;
        auto raw = SampleQuery(full[i]);
        if (raw.cloud->size() < 200 || !HasImuCoverage(raw,imu,maximum_shift+.001)) return;
        auto target = ReferenceMap(corrected[i-1],corrected[i+1],fits[i],fits[i+1]);
        DeskewState state;
        state.acceleration_scale=9.81/param.initialization_acceleration_norm;
        {
            auto next = std::lower_bound(states.begin(),states.end(),raw.time,
                [](const auto& s,double t){return s.time<t;});
            if (next==states.end() || (next==states.begin() && next->time>raw.time)) return;
            const auto& b=*next;
            const auto& a=next==states.begin() ? b : *(next-1);
            if (b.time-a.time>.2) return;
            const double u=b.time>a.time ? (raw.time-a.time)/(b.time-a.time) : 0.;
            const auto rotation=Eigen::Quaterniond(a.rotation).slerp(u,Eigen::Quaterniond(b.rotation)).toRotationMatrix();
            state.velocity=rotation.transpose()*((1-u)*a.velocity+u*b.velocity);
            state.gravity=rotation.transpose()*((1-u)*a.gravity+u*b.gravity);
            state.gyro_bias=(1-u)*a.gyro_bias+u*b.gyro_bias;
            state.acceleration_bias=(1-u)*a.acceleration_bias+u*b.acceleration_bias;
        }
        prepared_queries[job].emplace(Query{i,std::move(raw),FrozenPlaneMap(target),state,full[i].time >= split});
    });
    // Keep query admission and order identical to the original sequential loop.
    for (auto& query : prepared_queries) if (query) queries.push_back(std::move(*query));
    std::ofstream details(diagnostics/"imu_time_offset_matches.csv");
    std::ofstream curve(diagnostics/"imu_time_offset_curve.csv");
    if (!details || !curve) throw std::runtime_error("cannot create IMU offset matching diagnostics");
    details << "scan,relative_time,split,stamp_shift_ms,effective_imu_time_offset_ms,registered,points,matches,match_ratio,inlier_rms_m,capped_rms_m,balanced_score\n" << std::setprecision(17);
    curve << "stamp_shift_ms,effective_imu_time_offset_ms,train_match_ratio,validation_match_ratio,match_ratio,inlier_rms_m,capped_rms_m,points,matches,registered_scans,query_scans,train_balanced_score,validation_balanced_score,balanced_score\n" << std::setprecision(17);
    const auto scoring_start=std::chrono::steady_clock::now();
    std::vector<livo::OffsetCandidateScore> scores;
    std::vector<livo::OffsetMatchScore> query_scores(queries.size());
    for (double shift : candidate_shifts) {
        if(param.progress)param.progress();
        livo::OffsetCandidateScore candidate; candidate.stamp_shift_s = shift;
        ForEachInitializationTask(queries.size(),param.workers,[&](size_t i) {
            const auto& query=queries[i];
            const auto cloud=DeskewWithImu(query.raw,imu,query.state,candidate.stamp_shift_s);
            query_scores[i]=query.map.Evaluate(cloud);
        });
        // Callback, accumulation and CSV writes remain on the caller, in the
        // original query order. Worker completion order cannot change a score.
        for(size_t i=0;i<queries.size();++i) {
            if(param.progress)param.progress();
            const auto& query=queries[i];
            const auto& score=query_scores[i];
            (query.validation ? candidate.validation : candidate.train) += score;
            details << query.index << ',' << query.raw.time << ',' << (query.validation ? "validation" : "train")
                    << ',' << candidate.stamp_shift_s*1000 << ',' << (param.imu_time_offset_s-candidate.stamp_shift_s)*1000
                    << ',' << score.valid_scans << ',' << score.points << ',' << score.matches << ',' << score.ratio()
                    << ',' << score.rms() << ',' << score.capped_rms() << ',';
            if (score.valid()) details << Policy::BalancedScore(score);
            details << '\n';
        }
        const auto total = candidate.total();
        curve << candidate.stamp_shift_s*1000 << ',' << (param.imu_time_offset_s-candidate.stamp_shift_s)*1000
              << ',' << candidate.train.ratio() << ',' << candidate.validation.ratio() << ',' << total.ratio()
              << ',' << total.rms() << ',' << total.capped_rms() << ',' << total.points << ',' << total.matches
              << ',' << total.valid_scans << ',' << total.scans << ',';
        if (candidate.train.valid()) curve << Policy::BalancedScore(candidate.train);
        curve << ',';
        if (candidate.validation.valid()) curve << Policy::BalancedScore(candidate.validation);
        curve << ',';
        if (total.valid()) curve << Policy::BalancedScore(total);
        curve << '\n';
        scores.push_back(candidate);
        std::clog << "[imu time offset init] candidate stamp_shift_ms=" << candidate.stamp_shift_s*1000
                       << " match_ratio=" << total.ratio() << " rms_m=" << total.rms()
                       << " balanced_score=" << Policy::BalancedScore(total);
    }
    if (!details || !curve) throw std::runtime_error("cannot write IMU offset matching diagnostics");
    const auto scoring_end=std::chrono::steady_clock::now();
    Json::Value timing;
    timing["workers"]=param.workers;
    timing["reference_registration_and_deskew_s"]=std::chrono::duration<double>(maps_start-reference_start).count();
    timing["frozen_query_maps_s"]=std::chrono::duration<double>(scoring_start-maps_start).count();
    timing["candidate_scoring_s"]=std::chrono::duration<double>(scoring_end-scoring_start).count();
    std::ofstream(diagnostics/"imu_time_offset_timing.json")<<timing<<'\n';
    auto out = Policy::Estimate(std::move(scores),param.imu_time_offset_s);
    out.scans = full.size(); out.window_s = last-begin;
    out.window_begin_s = begin; out.window_end_s = last;
    if (out.applied) {
        // A numerical maximum alone is not evidence for a clock correction.
        const double gain=Policy::BalancedScore(out.after)-Policy::BalancedScore(out.before);
        const double heldout_gain=Policy::BalancedScore(out.validation_after)-Policy::BalancedScore(out.validation_before);
        if (out.window_s<2.7 || out.observations<5 || !std::isfinite(gain) || gain<.001 ||
            !std::isfinite(heldout_gain) || heldout_gain<0.) {
            out.applied=false; out.estimated=false; out.stamp_shift_s=0.;
            out.after=out.before; out.validation_after=out.validation_before;
            out.reason="insufficient_offset_evidence_keep_configured_clock";
        }
    }
    return out;
}

Json::Value ImuTimeOffsetJson(const livo::ImuTimeOffsetResult& r,
                             double configured_offset_s) {
    Json::Value out;
    out["requested"] = r.requested; out["offset_estimated"] = r.estimated;
    out["applied"] = r.applied; out["reason"] = r.reason;
    out["method"] = "imu_deskew_frozen_map_balanced_score_search";
    out["policy"] = "two_stage_bias_motion_offset_initialization";
    out["calibrated_sensor"] = "imu"; out["reference_sensor"] = "lidar";
    out["configured_imu_time_offset_s"] = configured_offset_s;
    out["estimated_stamp_shift_ms"] = r.stamp_shift_s*1000;
    out["candidate_stamp_shift_ms"] = r.candidate_shift_s*1000;
    out["effective_imu_time_offset_s"] = configured_offset_s-r.stamp_shift_s;
    out["effective_stamp_shift_ms"] = (-configured_offset_s+r.stamp_shift_s)*1000;
    out["sign_convention"] = "corrected_imu_time = raw_time - configured_offset + stamp_shift";
    out["window_s"] = r.window_s; out["scans"] = Json::UInt64(r.scans);
    out["requested_window_s"] = Policy::kWindowSeconds;
    out["window_begin_s"] = r.window_begin_s; out["window_end_s"] = r.window_end_s;
    out["state_source"] = "frozen_bias_velocity_gravity_from_stage_two_livo";
    out["query_scans"] = Json::UInt64(r.observations);
    out["validation_query_scans"] = Json::UInt64(r.validation_observations);
    out["validation_stamp_shift_ms"] = r.validation_shift_s*1000;
    out["before"] = ScoreJson(r.before); out["after"] = ScoreJson(r.after);
    out["validation_before"] = ScoreJson(r.validation_before);
    out["validation_after_selected_offset"] = ScoreJson(r.validation_after);
    out["reference_map"] = "two_neighbor_lidar_only_deskewed_scans_query_excluded_frozen_for_all_candidates";
    out["deskew"] = "midpoint_integration_with_frozen_bias_corrected_stage_two_state";
    out["scoring"] = "final_pose_reassociation_fixed_points_and_gates_maximum_balanced_score";
    out["balanced_score_formula"] = "match_ratio * (1 - inlier_rms_m^2 / match_distance_m^2)";
    out["physical_clock_ground_truth"] = false;
    out["search_min_imu_time_offset_ms"] = Policy::kMinimumOffset*1000;
    out["search_max_imu_time_offset_ms"] = Policy::kMaximumOffset*1000;
    out["grid_step_ms"] = Policy::kGridStep*1000;
    out["match_distance_m"] = Policy::kMatchDistance;
    out["residual_cap_m"] = Policy::kResidualCap;
    out["selection_policy"] = "highest_balanced_score_with_minimum_gain_and_holdout_non_degradation";
    out["temporal_subsets_policy"] = "holdout_must_not_degrade_for_clock_change";
    out["candidates"] = Json::arrayValue;
    for (const auto& c : r.candidates) {
        Json::Value candidate;
        candidate["stamp_shift_ms"] = c.stamp_shift_s*1000;
        candidate["effective_imu_time_offset_s"] = configured_offset_s-c.stamp_shift_s;
        candidate["in_search_range"] = Policy::InSearchRange(configured_offset_s-c.stamp_shift_s);
        candidate["train"] = ScoreJson(c.train); candidate["validation"] = ScoreJson(c.validation);
        candidate["total"] = ScoreJson(c.total()); out["candidates"].append(candidate);
    }
    return out;
}
}  // namespace initialization
