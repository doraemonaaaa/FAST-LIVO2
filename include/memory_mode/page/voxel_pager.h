// Adapted from deep_autolabel frontend/fast-livo2/page/voxel_pager.h.
#pragma once
#include <filesystem>
#include <fstream>
#include <memory>
#include <unordered_set>
#include "voxel_map.h"

namespace flivo {
// Ephemeral, same-process lossless storage, not a portable map file. Keep keys
// in the original hash table so eviction cannot change its iteration order.
class VoxelPager {
 public:
    using Map=std::unordered_map<VOXEL_LOCATION,VoxelOctoTree*>;
    using Keys=std::unordered_set<VOXEL_LOCATION>;
    struct Stats { size_t cold=0, working=0; uint64_t loads=0,evictions=0,bytes_written=0; } stats;
    VoxelPager(Map& map,const std::filesystem::path& path):map_(map),path_(path) {
        std::filesystem::create_directories(path.parent_path());
        file_.open(path,std::ios::binary|std::ios::in|std::ios::out|std::ios::trunc);
        if(!file_) throw std::runtime_error("cannot create voxel backing file");
        for(const auto& item:map_)if(item.second)resident_.insert(item.first);
    }
    ~VoxelPager() {file_.close();std::error_code ec;std::filesystem::remove(path_,ec);}
    void BeginFrame() { previous_working_.swap(frame_working_);frame_working_.clear(); }
    VoxelOctoTree* Ensure(const VOXEL_LOCATION& key,bool dirty=false,bool track=true) {
        auto it=map_.find(key);if(it==map_.end()) return nullptr;
        if(!it->second) {
            auto entry=entries_.find(key);
            if(entry==entries_.end() || !entry->second.saved) throw std::runtime_error("missing cold voxel state");
            file_.flush();file_.clear();file_.seekg(entry->second.offset);
            auto node=ReadNode(0);if(!file_) throw std::runtime_error("cannot restore cold voxel");
            it->second=node.release();--stats.cold;++stats.loads;
        }
        resident_.insert(key);if(track)frame_working_.insert(key);
        if(dirty) entries_[key].dirty=true;
        return it->second;
    }
    void Prepare(const std::vector<pointWithVar>& points,double size) {
        Keys needed;needed.reserve(points.size()*2);
        for(const auto& p:points) {
            const auto key=Location(p.point_w,size);int dir[3];const int64_t cells[3]={key.x,key.y,key.z};
            for(int k=0;k<3;++k) {const double u=p.point_w[k]/size-cells[k];dir[k]=u<.25?-1:(u>.75?1:0);}
            for(int x=0;x<=(dir[0]!=0);++x)for(int y=0;y<=(dir[1]!=0);++y)for(int z=0;z<=(dir[2]!=0);++z)
                needed.emplace(key.x+x*dir[0],key.y+y*dir[1],key.z+z*dir[2]);
        }
        stats.working=0;for(const auto& key:needed)stats.working+=map_.count(key);
        // Retain actual accesses from this and the preceding LIO frame, not a
        // guessed radius or point-count cap. Recompute on every EKF iteration.
        frame_working_.insert(needed.begin(),needed.end());
        needed.insert(frame_working_.begin(),frame_working_.end());
        needed.insert(previous_working_.begin(),previous_working_.end());
        for(auto it=resident_.begin();it!=resident_.end();) {
            if(!needed.count(*it)) {auto item=map_.find(*it);Evict(*item);it=resident_.erase(it);}
            else ++it;
        }
        for(const auto& key:needed) Ensure(key,false,false);
    }
    void Forget(const VOXEL_LOCATION& key) {
        auto it=map_.find(key);if(it!=map_.end() && !it->second)--stats.cold;
        entries_.erase(key);resident_.erase(key);frame_working_.erase(key);previous_working_.erase(key);
    }
    static VOXEL_LOCATION Location(const V3D& p,double s) {
        return {static_cast<int64_t>(std::floor(p.x()/s)),static_cast<int64_t>(std::floor(p.y()/s)),static_cast<int64_t>(std::floor(p.z()/s))};
    }
 private:
    struct Entry {std::streamoff offset=0;uint64_t capacity=0;bool saved=false,dirty=true;};
    template<class T> void Write(const T& value) {if(sizing_)record_size_+=sizeof(value);else file_.write(reinterpret_cast<const char*>(&value),sizeof(value));}
    template<class T> void Read(T& value) {file_.read(reinterpret_cast<char*>(&value),sizeof(value));}
    template<class M> void WriteMatrix(const M& value) {const auto bytes=value.size()*sizeof(typename M::Scalar);if(sizing_)record_size_+=bytes;else file_.write(reinterpret_cast<const char*>(value.data()),bytes);}
    template<class M> void ReadMatrix(M& value) {file_.read(reinterpret_cast<char*>(value.data()),value.size()*sizeof(typename M::Scalar));}
    void Point(pointWithVar& p,bool read) {
        auto matrix=[&](auto& v){if(read)ReadMatrix(v);else WriteMatrix(v);};
        matrix(p.point_b);matrix(p.point_i);matrix(p.point_w);matrix(p.var_nostate);matrix(p.body_var);matrix(p.var);matrix(p.point_crossmat);matrix(p.normal);
    }
    void Plane(VoxelPlane& p,bool read) {
        auto matrix=[&](auto& v){if(read)ReadMatrix(v);else WriteMatrix(v);};
        auto scalar=[&](auto& v){if(read)Read(v);else Write(v);};
        matrix(p.center_);matrix(p.normal_);matrix(p.y_normal_);matrix(p.x_normal_);matrix(p.covariance_);matrix(p.plane_var_);
        scalar(p.radius_);scalar(p.min_eigen_value_);scalar(p.mid_eigen_value_);scalar(p.max_eigen_value_);scalar(p.d_);scalar(p.points_size_);scalar(p.is_plane_);scalar(p.is_init_);scalar(p.id_);scalar(p.is_update_);
    }
    void Scalars(VoxelOctoTree& n,bool read) {
        auto s=[&](auto& v){if(read)Read(v);else Write(v);};
        s(n.layer_);s(n.octo_state_);for(auto& v:n.voxel_center_)s(v);s(n.quater_length_);s(n.planer_threshold_);s(n.points_size_threshold_);s(n.update_size_threshold_);s(n.max_points_num_);s(n.max_layer_);s(n.new_points_);s(n.init_octo_);s(n.update_enable_);
    }
    void WriteNode(VoxelOctoTree& n) {
        Scalars(n,false);Plane(*n.plane_ptr_,false);
        uint64_t count=n.layer_init_num_.size();Write(count);for(auto& v:n.layer_init_num_)Write(v);
        count=n.temp_points_.size();Write(count);for(auto& p:n.temp_points_)Point(p,false);
        for(auto* child:n.leaves_) {const bool present=child;Write(present);if(child)WriteNode(*child);}
    }
    std::unique_ptr<VoxelOctoTree> ReadNode(int depth) {
        if(depth>6) throw std::runtime_error("invalid voxel cache depth");
        auto n=std::make_unique<VoxelOctoTree>(0,0,5,100,.01f);Scalars(*n,true);Plane(*n->plane_ptr_,true);
        uint64_t count=0;Read(count);if(count>16)throw std::runtime_error("invalid voxel layer cache");
        n->layer_init_num_.resize(count);for(auto& v:n->layer_init_num_)Read(v);
        Read(count);if(count>1000000)throw std::runtime_error("invalid voxel point cache");
        n->temp_points_.resize(count);for(auto& p:n->temp_points_)Point(p,true);
        for(auto& child:n->leaves_) {bool present=false;Read(present);if(present)child=ReadNode(depth+1).release();}
        if(!file_)throw std::runtime_error("truncated voxel backing record");
        return n;
    }
    void Evict(Map::value_type& item) {
        auto& entry=entries_[item.first];
        if(!entry.saved || entry.dirty) {
            // Reuse a voxel's slot once its state size stabilizes; do not append
            // another full version of the same voxel on every frame.
            sizing_=true;record_size_=0;WriteNode(*item.second);sizing_=false;
            const bool reuse=entry.saved && record_size_<=entry.capacity;
            file_.clear();if(reuse)file_.seekp(entry.offset);else file_.seekp(0,std::ios::end);
            const auto begin=file_.tellp();WriteNode(*item.second);
            if(!reuse) {
                uint64_t capacity=4096;while(capacity<record_size_)capacity*=2;
                if(capacity>record_size_) {file_.seekp(begin+static_cast<std::streamoff>(capacity-1));file_.put(0);}
                entry.capacity=capacity;
            }
            if(!file_)throw std::runtime_error("cannot write cold voxel");
            entry.offset=begin;
            entry.saved=true;entry.dirty=false;stats.bytes_written+=record_size_;
        }
        delete item.second;item.second=nullptr;++stats.evictions;++stats.cold;
    }
    bool sizing_=false;uint64_t record_size_=0;
    Keys resident_,frame_working_,previous_working_;
    Map& map_;std::filesystem::path path_;std::fstream file_;std::unordered_map<VOXEL_LOCATION,Entry> entries_;
};
}
