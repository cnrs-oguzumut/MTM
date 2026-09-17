// AvalancheRecorder.h
#ifndef AVALANCHE_RECORDER_H
#define AVALANCHE_RECORDER_H

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <cmath>

enum class AvalancheEventType {
    INITIAL_GUESS,
    LBFGS_ITER,
    LBFGS_CONVERGED,
    BEFORE_REMESH,
    AFTER_REMESH,
    REMESH_REJECTED,
    REMESH_ACCEPTED,
    AVALANCHE_COMPLETE
};

inline std::string to_string(AvalancheEventType type) {
    switch (type) {
        case AvalancheEventType::INITIAL_GUESS: return "INITIAL_GUESS";
        case AvalancheEventType::LBFGS_ITER: return "LBFGS_ITER";
        case AvalancheEventType::LBFGS_CONVERGED: return "LBFGS_CONVERGED";
        case AvalancheEventType::BEFORE_REMESH: return "BEFORE_REMESH";
        case AvalancheEventType::AFTER_REMESH: return "AFTER_REMESH";
        case AvalancheEventType::REMESH_REJECTED: return "REMESH_REJECTED";
        case AvalancheEventType::REMESH_ACCEPTED: return "REMESH_ACCEPTED";
        case AvalancheEventType::AVALANCHE_COMPLETE: return "AVALANCHE_COMPLETE";
        default: return "UNKNOWN";
    }
}

struct MicroStepRecord {
    int global_micro_step = 0;
    int load_step = 0;
    std::string phase = "INITIAL_RELAX"; // "INITIAL_RELAX", "REMESH_PASS_1", etc.
    AvalancheEventType event_type = AvalancheEventType::LBFGS_ITER;
    int lbfgs_iter = -1;
    double energy = 0.0;
    double energy_change = 0.0;
    double stress = 0.0;
    double grad_norm = 0.0;
    int num_elements = 0;
    int triangles_changed = 0;
    std::string note = "";
};

class AvalancheRecorder {
private:
    bool enabled_ = false;
    int current_load_step_ = -1;
    std::string current_phase_ = "INITIAL_RELAX";
    int global_micro_counter_ = 0;
    double last_energy_ = 0.0;
    std::vector<MicroStepRecord> buffer_;

public:
    static AvalancheRecorder& instance() {
        static AvalancheRecorder recorder;
        return recorder;
    }

    void setEnabled(bool enabled) {
        enabled_ = enabled;
    }

    bool isEnabled() const {
        return enabled_;
    }

    void startLoadStep(int load_step) {
        if (!enabled_) return;
        current_load_step_ = load_step;
        current_phase_ = "INITIAL_RELAX";
        global_micro_counter_ = 0;
        last_energy_ = 0.0;
        buffer_.clear();
    }

    void setPhase(const std::string& phase) {
        current_phase_ = phase;
    }

    void recordStep(AvalancheEventType event_type,
                    int lbfgs_iter,
                    double energy,
                    double stress,
                    double grad_norm,
                    int num_elements = 0,
                    int triangles_changed = 0,
                    const std::string& note = "") {
        if (!enabled_) return;

        double de = buffer_.empty() ? 0.0 : (energy - last_energy_);
        last_energy_ = energy;

        MicroStepRecord rec;
        rec.global_micro_step = global_micro_counter_++;
        rec.load_step = current_load_step_;
        rec.phase = current_phase_;
        rec.event_type = event_type;
        rec.lbfgs_iter = lbfgs_iter;
        rec.energy = energy;
        rec.energy_change = de;
        rec.stress = stress;
        rec.grad_norm = grad_norm;
        rec.num_elements = num_elements;
        rec.triangles_changed = triangles_changed;
        rec.note = note;

        buffer_.push_back(rec);
    }

    bool hasRecords() const {
        return !buffer_.empty();
    }

    void commitToFile(const std::string& directory = "avalanche_trace",
                      int pre_file_id = -1,
                      int post_file_id = -1) {
        if (!enabled_ || buffer_.empty()) return;

        std::filesystem::create_directories(directory);

        std::stringstream ss;
        if (pre_file_id >= 0 && post_file_id >= 0) {
            ss << directory << "/config_" << std::setw(5) << std::setfill('0') << pre_file_id
               << "_to_" << std::setw(5) << std::setfill('0') << post_file_id
               << "_step_" << std::setw(5) << std::setfill('0') << current_load_step_ << ".csv";
        } else {
            ss << directory << "/step_" << std::setw(5) << std::setfill('0') << current_load_step_ << "_surgery.csv";
        }
        std::string filename = ss.str();

        std::ofstream file(filename);
        if (!file) {
            std::cerr << "AvalancheRecorder: could not open " << filename << " for writing" << std::endl;
            return;
        }

        file << "global_micro_step,load_step,phase,event_type,lbfgs_iter,energy,energy_change,stress,grad_norm,num_elements,triangles_changed,note\n";
        file << std::scientific << std::setprecision(12);

        for (const auto& rec : buffer_) {
            file << rec.global_micro_step << ","
                 << rec.load_step << ","
                 << rec.phase << ","
                 << to_string(rec.event_type) << ","
                 << rec.lbfgs_iter << ","
                 << rec.energy << ","
                 << rec.energy_change << ","
                 << rec.stress << ","
                 << rec.grad_norm << ","
                 << rec.num_elements << ","
                 << rec.triangles_changed << ","
                 << "\"" << rec.note << "\"\n";
        }
        file.close();
        std::cout << "✓ Avalanche surgery trace written: " << filename
                  << " (" << buffer_.size() << " micro-steps)" << std::endl;
    }

    void discard() {
        buffer_.clear();
    }
};

#endif // AVALANCHE_RECORDER_H
