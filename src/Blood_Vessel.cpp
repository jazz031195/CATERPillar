#include "Blood_Vessel.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>


using namespace Eigen;
using namespace std;


Blood_Vessel::Blood_Vessel(){}

Blood_Vessel::~Blood_Vessel()
{
    ramification_spheres.clear();
}


void Blood_Vessel::keep_one_sphere(){
    Sphere s (this->ramification_spheres[0][0]);
    ramification_spheres[0].clear();
    add_sphere(s);
    growth_attempts += 1;


}

void Blood_Vessel::destroy(){
    ramification_spheres.clear();
    lengths_branches.clear();
    attractors.clear();
    children_branches.clear();
    growth_attempts = 0;

}

void Blood_Vessel::addToGrid(SphereGrid &grid) const{

    for (const auto &branch : ramification_spheres){
        for (const auto &sph : branch){
            grid.insert(sph);
        }
    }
}

void Blood_Vessel::add_first_sphere(const Sphere &s){

    Sphere s0 = s;
    s0.branch_id = 0;
    ramification_spheres.clear();
    ramification_spheres.resize(1);
    ramification_spheres[0].push_back(s0);

}


void Blood_Vessel::add_sphere(const Sphere &sphere_to_add){
    // add sphere to the main vessel (branch 0)
    this->ramification_spheres[0].push_back(sphere_to_add);
}


void Blood_Vessel::update_Volume(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits){
    double new_volume = 0.0;

    const std::vector<Sphere> &trunk = ramification_spheres[0];

    for (size_t i = 1; i < trunk.size(); ++i) {
        const Sphere &last_sphere = trunk[i - 1];
        const Sphere &current_sphere = trunk[i];

        double distance = (current_sphere.center - last_sphere.center).norm();

        // Volume of truncated cone
        double segment_volume = M_PI * distance * (
            current_sphere.radius * current_sphere.radius +
            last_sphere.radius * last_sphere.radius +
            current_sphere.radius * last_sphere.radius) / 3.0;

        bool current_in_bounds = Obstacle::check_borders(min_limits, max_limits, current_sphere.center, barrier_tickness);
        bool last_in_bounds = Obstacle::check_borders(min_limits, max_limits, last_sphere.center, barrier_tickness);

        if (current_in_bounds && last_in_bounds) {
            new_volume += segment_volume;
        } else if (current_in_bounds || last_in_bounds) {
            new_volume += segment_volume / 2.0;
        }
        // If both are out of bounds, add nothing
    }

    volume = new_volume;
}

void Blood_Vessel::compute_processes_icvf(const int &factor, const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits) {
    volume_processes = 0.0;
    for (int b = 1; b < ramification_spheres.size(); ++b) {
        for (int i = factor; i < ramification_spheres[b].size(); i += factor) {
            Sphere &s = ramification_spheres[b][i-factor];
            Sphere &s_next = ramification_spheres[b][i];
            if (s_next.center [0] + s_next.radius < min_limits[0] || s_next.center [0] - s_next.radius > max_limits[0] ||
                s_next.center [1] + s_next.radius < min_limits[1] || s_next.center [1] - s_next.radius > max_limits[1] ||
                s_next.center [2] + s_next.radius < min_limits[2] || s_next.center [2] - s_next.radius > max_limits[2]) {
                continue;
            }

            double distance = (s_next.center - s.center).norm();
            double v = M_PI * (s.radius * s.radius + s_next.radius * s_next.radius + s.radius * s_next.radius) * distance / 3.0;
            volume_processes += v;

        }
    }
}

bool Blood_Vessel::findSphereById(const int &sphere_id, int &out_branch, int &out_index) const {
    for (int b = 0; b < ramification_spheres.size(); ++b) {
        for (int i = 0; i < ramification_spheres[b].size(); ++i) {
            if (ramification_spheres[b][i].id == sphere_id) {
                out_branch = b;
                out_index = i;
                return true;
            }
        }
    }
    return false;
}

void Blood_Vessel::enforceMurraysLaw() {
    // Walk the branch tree from the trunk (branch 0) outward via children_branches,
    // so every parent branch is fully corrected before any of its children are
    // processed (a correction at one junction must cascade to junctions further
    // downstream on the same branch).
    std::vector<int> queue = {0};
    for (size_t qi = 0; qi < queue.size(); ++qi) {
        int parent_branch = queue[qi];
        if (parent_branch >= static_cast<int>(children_branches.size())) {
            continue;
        }

        // Find where each child attaches along the parent, and process them in
        // that order (most upstream first). A correction only rescales the parent
        // *downstream* of its junction, so processing upstream-to-downstream means
        // an already-finalized (more upstream) sibling junction is never touched
        // again by a later (more downstream) sibling's correction.
        std::vector<std::pair<int, int>> children_by_position; // (parent_index, child_branch)
        for (int child_branch : children_branches[parent_branch]) {
            queue.push_back(child_branch);

            if (child_branch >= static_cast<int>(ramification_spheres.size()) || ramification_spheres[child_branch].empty()) {
                continue;
            }

            int parent_index = -1;
            int child_parent_id = ramification_spheres[child_branch][0].parent_id;
            for (int i = 0; i < static_cast<int>(ramification_spheres[parent_branch].size()); ++i) {
                if (ramification_spheres[parent_branch][i].id == child_parent_id) {
                    parent_index = i;
                    break;
                }
            }
            if (parent_index < 0) {
                continue; // no matching parent sphere found; leave this branch untouched
            }
            children_by_position.push_back({parent_index, child_branch});
        }
        std::sort(children_by_position.begin(), children_by_position.end());

        for (const auto &entry : children_by_position) {
            int parent_index = entry.first;
            int child_branch = entry.second;

            double R_parent = ramification_spheres[parent_branch][parent_index].radius;

            // current radii of the two daughters just downstream of the junction: the
            // branch that just started (r2), and the parent's own continuation (r1),
            // which by default both start out equal to R_parent.
            double r2 = ramification_spheres[child_branch][0].radius;
            bool has_continuation = parent_index + 1 < static_cast<int>(ramification_spheres[parent_branch].size());
            double r1 = has_continuation ? ramification_spheres[parent_branch][parent_index + 1].radius : R_parent;

            double sum_cubes = r1 * r1 * r1 + r2 * r2 * r2;
            if (sum_cubes <= 0.0) {
                continue;
            }

            // Rescale both daughters by the same factor so R_parent^3 = r1^3 + r2^3 holds;
            // since r1 and r2 both start out equal to R_parent, this splits evenly.
            double scale = std::cbrt((R_parent * R_parent * R_parent) / sum_cubes);

            // Predict whether the child would survive this scale *before* applying
            // it anywhere. A branch's radius decays monotonically from its start to
            // its tip, so the tip is always the weakest link: if the scaled tip
            // radius still clears minimum_radius, every other sphere in the branch
            // does too, and the branch is guaranteed to still reach the voxel plane
            // it was grown to cross (scaling never moves sphere centers). If the
            // scaled tip would drop below minimum_radius, the branch would end up
            // truncated downstream anyway -- almost always losing exactly the part
            // that crossed the plane -- so discard it outright here and skip the
            // split entirely, leaving the parent's own radius untouched, as though
            // this branch had never been created. This avoids needlessly thinning
            // a parent for a branch that ultimately doesn't survive.
            double predicted_child_tip_radius = ramification_spheres[child_branch].back().radius * scale;

            // Symmetrically, check the parent's own continuation. Unlike a branch,
            // the trunk's radius isn't a strictly monotonic decay (it can "bead" --
            // oscillate around a baseline), so its weakest point downstream of this
            // junction isn't necessarily its very last sphere: scan the whole
            // affected range for the true minimum. A parent that has already
            // received several upstream splits can be pushed below minimum_radius
            // here even when this specific split looks reasonable in isolation --
            // in which case this split would only get truncated downstream anyway,
            // so reject it the same way: discard the child, leave the parent as-is.
            double predicted_parent_tail_radius = R_parent;
            if (has_continuation) {
                double min_parent_tail_radius = std::numeric_limits<double>::max();
                for (int i = parent_index + 1; i < static_cast<int>(ramification_spheres[parent_branch].size()); ++i) {
                    min_parent_tail_radius = std::min(min_parent_tail_radius, ramification_spheres[parent_branch][i].radius);
                }
                predicted_parent_tail_radius = min_parent_tail_radius * scale;
            }

            if (predicted_child_tip_radius < minimum_radius || predicted_parent_tail_radius < minimum_radius) {
                deleteSubtree(child_branch);
                continue;
            }

            if (has_continuation) {
                for (int i = parent_index + 1; i < static_cast<int>(ramification_spheres[parent_branch].size()); ++i) {
                    ramification_spheres[parent_branch][i].radius *= scale;
                }
            }
            for (auto &sph : ramification_spheres[child_branch]) {
                sph.radius *= scale;
            }
        }
    }
}

void Blood_Vessel::deleteSubtree(int branch) {
    if (branch < 0 || branch >= static_cast<int>(ramification_spheres.size())) {
        return;
    }
    ramification_spheres[branch].clear();
    if (branch < static_cast<int>(lengths_branches.size())) {
        lengths_branches[branch].clear();
    }
    if (branch < static_cast<int>(children_branches.size())) {
        for (int child : children_branches[branch]) {
            deleteSubtree(child);
        }
        children_branches[branch].clear();
    }
}

void Blood_Vessel::pruneUndersizedBranches(const Eigen::Vector3d &min_limits, const Eigen::Vector3d &max_limits) {
    // Walk the branch tree from the trunk (branch 0) outward. enforceMurraysLaw can
    // shrink a branch's own starting radius, or a parent's continuation partway
    // along its length -- either way, once a sphere drops below minimum_radius,
    // truncate that branch right there (everything from that sphere onward,
    // including any children rooted on the removed tail, is deleted with it). A
    // non-trunk branch left with fewer than 2 spheres is pointless and is deleted
    // outright.
    std::vector<int> queue = {0};
    for (size_t qi = 0; qi < queue.size(); ++qi) {
        int b = queue[qi];
        if (b >= static_cast<int>(ramification_spheres.size()) || ramification_spheres[b].empty()) {
            continue;
        }

        auto &spheres = ramification_spheres[b];
        int cutoff = -1;
        for (int i = 0; i < static_cast<int>(spheres.size()); ++i) {
            if (spheres[i].radius < minimum_radius) {
                cutoff = i;
                break;
            }
        }

        if (cutoff >= 0) {
            // delete any children attached at or beyond the sphere being cut off
            if (b < static_cast<int>(children_branches.size())) {
                std::vector<int> kept_children;
                for (int child : children_branches[b]) {
                    int parent_index = -1;
                    if (child < static_cast<int>(ramification_spheres.size()) && !ramification_spheres[child].empty()) {
                        int pid = ramification_spheres[child][0].parent_id;
                        for (int i = 0; i < static_cast<int>(spheres.size()); ++i) {
                            if (spheres[i].id == pid) { parent_index = i; break; }
                        }
                    }
                    if (parent_index < 0 || parent_index >= cutoff) {
                        deleteSubtree(child);
                    } else {
                        kept_children.push_back(child);
                    }
                }
                children_branches[b] = std::move(kept_children);
            }

            spheres.resize(cutoff);
            if (b < static_cast<int>(lengths_branches.size()) && lengths_branches[b].size() > spheres.size()) {
                lengths_branches[b].resize(spheres.size());
            }

            if (b != 0 && spheres.size() < 2) {
                deleteSubtree(b);
                continue; // nothing left of this branch to descend into
            }
        }

        // Branches only exist to cross a voxel plane: if truncation above cut off
        // the part that was outside the box (or, defensively, if it otherwise
        // no longer reaches), the surviving stub dead-ends inside the domain and
        // is discarded entirely rather than kept as-is.
        if (b != 0 && !spheres.empty()) {
            const Sphere &tip = spheres.back();
            bool reaches_plane = !Obstacle::check_borders(min_limits, max_limits, tip.center, tip.radius);
            if (!reaches_plane) {
                deleteSubtree(b);
                continue;
            }
        }

        if (b >= static_cast<int>(children_branches.size())) {
            continue;
        }
        for (int child : children_branches[b]) {
            queue.push_back(child);
        }
    }
}

void Blood_Vessel::bridgeJunctionGaps(const int &factor) {
    int next_id = 0;
    for (const auto &branch : ramification_spheres) {
        for (const auto &s : branch) {
            next_id = std::max(next_id, s.id + 1);
        }
    }

    for (size_t b = 1; b < ramification_spheres.size(); ++b) {
        if (ramification_spheres[b].empty()) continue;

        Sphere &child_first = ramification_spheres[b][0];
        int pb, pi;
        if (!findSphereById(child_first.parent_id, pb, pi)) continue;
        const Sphere &parent_sph = ramification_spheres[pb][pi];

        double dist = (child_first.center - parent_sph.center).norm();
        double r_min = std::min(child_first.radius, parent_sph.radius);
        double max_spacing = r_min / static_cast<double>(factor);
        if (max_spacing <= 0.0 || dist <= max_spacing) continue;

        bool has_lengths = b < lengths_branches.size() && lengths_branches[b].size() == ramification_spheres[b].size();
        double parent_length = (pb < static_cast<int>(lengths_branches.size()) && pi < static_cast<int>(lengths_branches[pb].size()))
                                    ? lengths_branches[pb][pi] : 0.0;
        double child_length = has_lengths ? lengths_branches[b][0] : 0.0;

        const int max_extra_spheres = 200;
        int nbr_extra = std::min(static_cast<int>(std::ceil(dist / max_spacing)) - 1, max_extra_spheres);

        std::vector<Sphere> bridge;
        std::vector<double> bridge_lengths;
        int prev_id = parent_sph.id;
        for (int k = 1; k <= nbr_extra; ++k) {
            double t = static_cast<double>(k) / (nbr_extra + 1);
            Eigen::Vector3d pos = parent_sph.center + t * (child_first.center - parent_sph.center);
            double rad = parent_sph.radius + t * (child_first.radius - parent_sph.radius);
            Sphere s(next_id++, child_first.object_id, child_first.object_type, pos, rad, static_cast<int>(b), prev_id);
            bridge.push_back(s);
            prev_id = s.id;
            if (has_lengths) {
                bridge_lengths.push_back(parent_length + t * (child_length - parent_length));
            }
        }

        child_first.parent_id = prev_id;
        ramification_spheres[b].insert(ramification_spheres[b].begin(), bridge.begin(), bridge.end());
        if (has_lengths) {
            lengths_branches[b].insert(lengths_branches[b].begin(), bridge_lengths.begin(), bridge_lengths.end());
        }
    }
}

void Blood_Vessel::reinterpolateAfterShrink(const int &factor) {
    // Find the highest sphere id used anywhere in this vessel, so newly inserted
    // spheres get ids that don't collide with any existing one (parent_id links
    // and children_branches lookups rely on ids being unique within the vessel).
    int next_id = 0;
    for (const auto &branch : ramification_spheres) {
        for (const auto &s : branch) {
            next_id = std::max(next_id, s.id + 1);
        }
    }

    for (int b = 0; b < static_cast<int>(ramification_spheres.size()); ++b) {
        auto &branch = ramification_spheres[b];
        if (branch.size() < 2) {
            continue;
        }

        bool has_lengths = b < static_cast<int>(lengths_branches.size()) && lengths_branches[b].size() == branch.size();

        std::vector<Sphere> new_branch;
        std::vector<double> new_lengths;
        new_branch.reserve(branch.size());
        if (has_lengths) new_lengths.reserve(branch.size());

        for (size_t i = 0; i + 1 < branch.size(); ++i) {
            new_branch.push_back(branch[i]);
            if (has_lengths) new_lengths.push_back(lengths_branches[b][i]);

            const Sphere &a = branch[i];
            double dist = (branch[i + 1].center - a.center).norm();
            // Size using the *smaller* of the two radii: radius interpolates linearly
            // between them, so the narrower end is where the r_max/factor constraint
            // is tightest. Sizing off it keeps every resulting local pair (checked
            // against its own local r_max) within bounds, not just the pair's original
            // endpoints.
            double r_min = std::min(a.radius, branch[i + 1].radius);
            double max_spacing = r_min / static_cast<double>(factor);

            if (max_spacing > 0.0 && dist > max_spacing) {
                // Safety cap: pruneUndersizedBranches should already keep radii from
                // collapsing toward zero, but a parent's own continuation isn't
                // pruned, so guard against a tiny max_spacing still making
                // dist/max_spacing astronomically large.
                const int max_extra_spheres_per_gap = 200;
                int nbr_extra = std::min(static_cast<int>(std::ceil(dist / max_spacing)) - 1, max_extra_spheres_per_gap);
                int prev_id = a.id;
                for (int k = 1; k <= nbr_extra; ++k) {
                    double t = static_cast<double>(k) / (nbr_extra + 1);
                    Eigen::Vector3d pos = a.center + t * (branch[i + 1].center - a.center);
                    double rad = a.radius + t * (branch[i + 1].radius - a.radius);
                    Sphere s(next_id++, a.object_id, a.object_type, pos, rad, a.branch_id, prev_id);
                    new_branch.push_back(s);
                    prev_id = s.id;
                    if (has_lengths) {
                        double len = lengths_branches[b][i] + t * (lengths_branches[b][i + 1] - lengths_branches[b][i]);
                        new_lengths.push_back(len);
                    }
                }
                // keep the parent chain consistent: the sphere after the gap now
                // attaches to the last sphere inserted, not directly to a.
                branch[i + 1].parent_id = prev_id;
            }
        }
        new_branch.push_back(branch.back());
        if (has_lengths) new_lengths.push_back(lengths_branches[b].back());

        branch = std::move(new_branch);
        if (has_lengths) {
            lengths_branches[b] = std::move(new_lengths);
        }
    }
}
