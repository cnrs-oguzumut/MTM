# Factor main.cpp into focused modules without changing calculations

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` from the repository root. That file requires this document to be self-contained, novice-readable, updated during implementation, and validated by observable behavior. The current working branch for this plan is `feature/main_factor`.

## Purpose / Big Picture

The file `src/main.cpp` currently mixes the program entry point, simulation examples, remeshing helpers, data analysis routines, acoustic tensor studies, indentation experiments, and utility functions in one large translation unit. This makes it hard to find the active experiment, risky to edit one experiment without touching another, and difficult to reuse helpers from new experiments.

After this refactor, a user should still be able to build and run the executable exactly as before, and the active `run_final_shift_tests(...)` workflow should produce the same files under `final_tests/`. The visible improvement is organizational: `src/main.cpp` should become a small entry-point file, while each family of experiments lives in a named `.cpp` file with a matching header. The calculations, function bodies, command-line arguments, output folder names, and saved data formats must not change.

## Progress

- [x] (2026-05-23T10:23Z) Created and switched to branch `feature/main_factor`.
- [x] (2026-05-23T10:23Z) Read `PLANS.md` and confirmed the required ExecPlan sections and formatting rules.
- [x] (2026-05-23T10:23Z) Mapped the current `src/main.cpp` top-level functions and build layout.
- [x] (2026-05-23T10:23Z) Authored this initial self-contained execution plan.
- [ ] Reconfigure and build the current branch before refactoring to establish a baseline.
- [ ] Create shared headers for functions that are currently cross-called inside `src/main.cpp`.
- [ ] Move utility and remeshing helper functions out of `src/main.cpp` without changing their bodies.
- [ ] Move each experiment family into focused source files, validating after each move.
- [ ] Reduce `src/main.cpp` to command-line parsing and dispatch.
- [ ] Run final validation and record the observed outputs in this plan.

## Surprises & Discoveries

- Observation: `src/main.cpp` still has 6,257 lines even after the recent `shift_vertical_horizontal` extraction.
  Evidence: `wc -l src/main.cpp src/shift_vertical_horizontal.cpp src/main_splines.cpp CMakeLists.txt` reported `6257 src/main.cpp`.

- Observation: `src/main.cpp` contains `#include` directives inside the file body, not only at the top. These are legal C++ preprocessor directives, but they hide dependencies and should be removed or moved to the correct new module include blocks during extraction.
  Evidence: `rg -n "^#include" src/main.cpp` showed includes at lines around 228, 3123, 3234, and 4963 as well as at the top of the file.

- Observation: The build uses `file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")`, filters out `main.cpp`, and explicitly adds `src/shift_vertical_horizontal.cpp`. New source files should be added explicitly or the CMake configure step must be rerun after adding them.
  Evidence: `CMakeLists.txt` around the `add_library(lattice_lib STATIC ...)` block.

- Observation: The working tree contains many untracked result folders and scripts. They are not part of this refactor and must not be deleted, reformatted, or committed unless the user explicitly asks.
  Evidence: `git status --short` lists untracked `.cache/`, `deneme*/`, `plots/`, `spectra/`, scripts, and `src/main_splines.cpp`.

## Decision Log

- Decision: Perform this as a pure refactor: move declarations and function bodies, but do not rewrite algorithms, change constants, rename output files, or alter data formats.
  Rationale: The user explicitly requested dividing `main.cpp` without losing features or changing calculations. Behavioral preservation is more important than making the code elegant in one pass.
  Date/Author: 2026-05-23 / Codex

- Decision: Use incremental module extraction with a build after each milestone, instead of one large move.
  Rationale: A large mechanical move can easily break hidden dependencies in this file. Incremental builds make it easier to identify the exact extraction that introduced a compile or link error.
  Date/Author: 2026-05-23 / Codex

- Decision: Keep public function names stable during the first refactor.
  Rationale: Existing scripts, local notes, or future manual calls may refer to names such as `example_1_conti_zanzotto`, `parametricAcousticStudy`, or `example_3_stress_controlled_final_clean`. Renaming can be a later cleanup with its own validation.
  Date/Author: 2026-05-23 / Codex

- Decision: Put experiment declarations under `include/experiments/` and source files under `src/experiments/` or `src/` depending on existing CMake simplicity.
  Rationale: The repository already has domain headers under `include/geometry`, `include/mesh`, `include/output`, and one new experiment header at `include/experiments/shift_vertical_horizontal.h`. Following that pattern makes the organization discoverable.
  Date/Author: 2026-05-23 / Codex

## Outcomes & Retrospective

No implementation milestone has been completed yet. The current outcome is this plan on branch `feature/main_factor`, ready for stepwise execution. This section must be updated after each major extraction and again after final validation.

## Context and Orientation

The repository root is `/Users/usalman/programming/FEM_2D/factorized/lattice_triangulation`. The project builds a C++ executable named `lattice_triangulation`. The active program entry point is `int main(int argc, char **argv)` in `src/main.cpp`. It currently parses optional command-line arguments `nx`, `ny`, `horizontal_steps`, and `vertical_steps`, then calls `run_final_shift_tests(...)`, declared in `include/experiments/shift_vertical_horizontal.h` and implemented in `src/shift_vertical_horizontal.cpp`.

A "translation unit" means one `.cpp` file after the C++ preprocessor has expanded its included headers. A "header" means a `.h` file that declares functions, classes, and types so other `.cpp` files can call them. This refactor should move implementation code from one large translation unit into several smaller translation units while preserving the same external behavior.

The important existing files are:

- `src/main.cpp`: currently contains many experiment functions and the executable `main`. This is the file to shrink.
- `src/shift_vertical_horizontal.cpp`: already contains the recently extracted final shift tests for the perturbed left-bottom cases.
- `include/experiments/shift_vertical_horizontal.h`: declares `run_final_shift_tests(...)`.
- `src/Remesher.cpp` and `include/mesh/Remesher.h`: contain the triangulation perturbation support used by the final shift tests.
- `CMakeLists.txt`: controls which `.cpp` files are compiled into the `lattice_lib` static library and the `lattice_triangulation` executable.
- `include/optimization/LatticeOptimizer.h` and `src/optimization/LatticeOptimizer.cpp`: define `UserData`, `map_points_to_solver_array`, `map_solver_array_to_points`, and `minimize_energy_with_triangles`, which many experiments call.
- `include/output/configuration_saver.h`, `include/defects/Defectanalysis.h`, `include/geometry/NeighborAnalyzer.h`, and the mesh/lattice headers: provide the simulation infrastructure used by the experiment functions.

At the time this plan was written, `src/main.cpp` top-level functions detected by search were:

- `example_1_atomistic_square()` inside a long commented block near the top.
- `debug_deformation_tests()`.
- `debug_deformation_tests_triangular()`.
- `example_2_conti_zanzotto_triangular()`.
- `hasConnectivityChanged(...)`.
- `writeSizesToFile(...)`.
- `memory(...)`.
- `example_1_conti_zanzotto(...)`.
- `analyze_data_from_folder(...)`.
- `example_1_shifting(...)`.
- `sortByCoordinates(...)`.
- `callback(...)`.
- `savePositionEtaStress(...)`.
- `single_dislo_LJ()`.
- `indentation()`.
- `parametricAcousticStudy()`.
- `parametricAcousticStudy_v2()`.
- `example_3_stress_controlled_final_clean(...)`.
- `main(...)`.

Some functions may have return types split across multiple lines and may not appear in the simple search result. The implementer must search before each move and not assume the above list is exhaustive.

## Plan of Work

Start by establishing a baseline on branch `feature/main_factor`. Run CMake configure and build with the current code. Do not run full simulations as the baseline unless the user asks, because full remeshing cases can be long. The baseline acceptance is that the project compiles and links the `lattice_triangulation` executable.

Next, create a small shared header for the helper functions in `src/main.cpp` that are called by more than one experiment module. The first candidate is `include/experiments/common_simulation_helpers.h`. It should declare only stable helpers that need to cross source-file boundaries, such as `writeSizesToFile(...)`, `calculateShapeDerivatives(...)`, and `perform_remeshing_loop_reduction(...)`. Do not move a helper into a header as an inline implementation unless it is already header-only or trivially templated. Keep implementation in `.cpp` files.

Then extract low-level helper implementations into `src/experiments/common_simulation_helpers.cpp` or `src/simulation_helpers.cpp`. Move `hasConnectivityChanged(...)`, `writeSizesToFile(...)`, and any small helper that is not tied to one experiment family. Preserve each function body exactly except for include fixes required by the new file. After the move, include the new header wherever needed and build.

After helpers compile, move restart and Zanzotto-family examples into a focused module. Create `include/experiments/zanzotto_examples.h` and `src/experiments/zanzotto_examples.cpp`. Move `memory(...)`, `example_1_conti_zanzotto(...)`, and `example_2_conti_zanzotto_triangular(...)` there. If those functions need helper declarations, include `common_simulation_helpers.h`. Build immediately after this move.

Next, move data analysis into `include/experiments/data_analysis.h` and `src/experiments/data_analysis.cpp`. Move `analyze_data_from_folder(...)` and any local-only helpers that it exclusively uses. If `analyze_data_from_folder(...)` uses functions from other modules, prefer declaring those dependencies in the smallest relevant header rather than including unrelated experiment headers.

Then move deformation and shifting examples into `include/experiments/deformation_examples.h` and `src/experiments/deformation_examples.cpp`. Move `debug_deformation_tests()`, `debug_deformation_tests_triangular()`, and `example_1_shifting(...)` there. Do not merge this with the already extracted `shift_vertical_horizontal` module because the final shift tests are active and should remain isolated.

Then move dislocation and indentation experiments into `include/experiments/dislocation_indentation.h` and `src/experiments/dislocation_indentation.cpp`. Move `sortByCoordinates(...)`, `callback(...)`, `savePositionEtaStress(...)`, `single_dislo_LJ()`, and `indentation()`. The word "callback" is vague, so preserve the current name during this refactor and add a short comment in the new file explaining which optimizer or workflow calls it, if that is clear from the surrounding code.

Then move acoustic tensor studies into `include/experiments/acoustic_studies.h` and `src/experiments/acoustic_studies.cpp`. Move `parametricAcousticStudy()` and `parametricAcousticStudy_v2()`. These functions likely depend on `include/acoustic_tensor.h`, `include/FEMHessianAssembler.h`, and several lattice/mesh headers, so keep their includes local to the new `.cpp`.

Then move stress-controlled examples into `include/experiments/stress_controlled_examples.h` and `src/experiments/stress_controlled_examples.cpp`. Move `example_3_stress_controlled_final_clean(...)` there.

Finally, reduce `src/main.cpp` to normal includes, command-line parsing, and dispatch to `run_final_shift_tests(...)`. Remove duplicate in-body `#include` directives from `main.cpp` as their corresponding code leaves the file. Keep old inactive example calls as comments only if they are still useful, but prefer comments that reference the new header and function location.

For every new `.cpp` file, update `CMakeLists.txt` explicitly. The project currently uses a source glob, but explicit listing prevents missing-file surprises in existing configured build directories. The pattern already used for `src/shift_vertical_horizontal.cpp` can be extended.

## Concrete Steps

All commands in this section should be run from the repository root:

    cd /Users/usalman/programming/FEM_2D/factorized/lattice_triangulation

Confirm the branch:

    git branch --show-current

Expected output:

    feature/main_factor

Establish a clean build baseline. If the build directory does not exist, CMake will create it:

    cmake -S . -B /private/tmp/lattice_triangulation_build_check
    cmake --build /private/tmp/lattice_triangulation_build_check --target lattice_triangulation

Expected successful ending:

    [100%] Built target lattice_triangulation

Before each extraction, inspect the function and its local dependencies. For example:

    rg -n "void memory|void example_1_conti_zanzotto|void example_2_conti_zanzotto_triangular" src/main.cpp
    sed -n '1230,2215p' src/main.cpp

When creating files, use `apply_patch` or an editor. Do not use destructive commands to delete source until the build has passed with the moved copy. A safe extraction loop is:

1. Add the new header and `.cpp` file.
2. Copy the target functions into the new `.cpp`.
3. Include the new header in `src/main.cpp` if `main.cpp` still calls those functions.
4. Build and fix missing includes or declarations.
5. Remove the original function bodies from `src/main.cpp`.
6. Build again.
7. Update this ExecPlan `Progress`, `Surprises & Discoveries`, and `Decision Log` if anything unexpected happened.

After each source-file addition, reconfigure once so CMake sees the new file:

    cmake -S . -B /private/tmp/lattice_triangulation_build_check
    cmake --build /private/tmp/lattice_triangulation_build_check --target lattice_triangulation

Use `git diff --stat` and `git diff --check` after each milestone:

    git diff --stat
    git diff --check

Expected output from `git diff --check` is no output. Any trailing whitespace or conflict marker must be fixed before continuing.

At the end, inspect that `main.cpp` is small and only dispatches:

    wc -l src/main.cpp
    rg -n "run_final_shift_tests|int main" src/main.cpp

The exact final line count is not fixed, but it should be dramatically smaller than 6,257 lines. A good target is under 250 lines if all inactive experiment bodies have moved out.

## Validation and Acceptance

Validation must prove both compilation and behavioral preservation.

First, the project must build:

    cmake -S . -B /private/tmp/lattice_triangulation_build_check
    cmake --build /private/tmp/lattice_triangulation_build_check --target lattice_triangulation

Accept only a build ending with:

    [100%] Built target lattice_triangulation

Second, the active executable dispatch must remain unchanged. Running the executable should still call `run_final_shift_tests(nx, ny, horizontal_steps, vertical_steps)`. For a small smoke test, use a tiny lattice and minimal step counts in a temporary directory so output does not pollute the repository:

    mkdir -p /private/tmp/main_factor_smoke
    cd /private/tmp/main_factor_smoke
    /private/tmp/lattice_triangulation_build_check/lattice_triangulation 4 4 1 1

Expected early output should include both final test folders:

    === Running final test: left_bottom_perturbed_no_remesh_amp2 ===
    shift_vertical_horizontal left_bottom_positive boundary setup:
    initial triangulation perturbation: enabled
    === Running final test: left_bottom_perturbed_remesh_amp2 ===

The smoke test should create:

    /private/tmp/main_factor_smoke/final_tests/left_bottom_perturbed_no_remesh_amp2
    /private/tmp/main_factor_smoke/final_tests/left_bottom_perturbed_remesh_amp2

If the smoke test takes too long in the remeshing case, stop it with `pkill -x lattice_triangulation` only after confirming that startup and the first case dispatch are correct. Record that interruption in `Outcomes & Retrospective`.

Third, compare the active-case source behavior before and after the refactor when possible. A practical method is to save the stdout of a tiny no-remesh run before moving code and compare it after the move. Because this repository currently runs case 5 and 6 together, the implementer may temporarily set the final test list to only the no-remesh perturbed case for comparison, but must revert that temporary edit before committing. If any comparison uses a temporary edit, record it in this plan.

Acceptance criteria:

- `src/main.cpp` contains the entry point and dispatch only, not thousands of lines of experiment implementations.
- All moved functions preserve their names, signatures, and bodies except for required include and namespace adjustments.
- `run_final_shift_tests(...)` still runs the perturbed left-bottom no-remesh and remesh cases with amplitude 2.
- The triangulation perturbation remains `-1e-7` for the perturbed initial mesh.
- Build succeeds after a clean CMake configure.
- No untracked result folders or unrelated scripts are committed as part of this refactor.

## Idempotence and Recovery

The refactor should be safe to repeat in small steps. CMake configure and build commands are idempotent: running them again should not alter source files. Moving a function should be done in two phases, copy then delete, so a failed build can be recovered by restoring the original function body from `git diff` or `git restore -p src/main.cpp`.

Do not use `git reset --hard` or delete output folders. If a move goes wrong, use:

    git diff
    git restore -p <file>

to selectively undo only the incorrect hunks. If a new module causes link errors, check whether the function declaration in the header exactly matches the implementation and whether the new `.cpp` file is part of `lattice_lib` in `CMakeLists.txt`.

If a source file is added but the linker cannot find a symbol, rerun:

    cmake -S . -B /private/tmp/lattice_triangulation_build_check

before debugging further, because the build directory may not have been regenerated after the new file was added.

## Artifacts and Notes

Initial research transcript:

    $ wc -l src/main.cpp src/shift_vertical_horizontal.cpp src/main_splines.cpp CMakeLists.txt
        6257 src/main.cpp
         536 src/shift_vertical_horizontal.cpp
        5432 src/main_splines.cpp
         542 CMakeLists.txt
       12767 total

Top-level function search in `src/main.cpp`:

    $ rg -n "^(void|int|double|bool|Eigen::|std::tuple|std::pair|std::vector|alglib::|DomainInfo|MatrixXd|VectorXd) [A-Za-z_][A-Za-z0-9_]*\\(" src/main.cpp
    232:void debug_deformation_tests() {
    427:void debug_deformation_tests_triangular() {
    614:void example_2_conti_zanzotto_triangular() {
    1023:bool hasConnectivityChanged(...)
    1226:void writeSizesToFile(int Nx, int Ny) {
    1239:void memory(int caller_id, int nx, int ny, int restart_iteration) {
    1680:void example_1_conti_zanzotto(int caller_id, int nx, int ny) {
    2214:void analyze_data_from_folder(...)
    2651:void example_1_shifting(int caller_id, int nx, int ny) {
    3283:void sortByCoordinates(...)
    3457:void callback(...)
    3552:void savePositionEtaStress(...)
    3740:void single_dislo_LJ() {
    4291:void indentation() {
    4726:void parametricAcousticStudy() {
    4989:void parametricAcousticStudy_v2() {
    5649:void example_3_stress_controlled_final_clean(...) {
    6221:int main(int argc, char **argv) {

Current untracked files to avoid committing unless explicitly requested:

    .cache/
    Makefile
    PLANS.md
    cauchy/
    deneme/
    deneme2/
    deneme3/
    dos_calculation.py
    logbin.py
    plot_stability.py
    plots/
    profile_output.txt
    profiling.py
    spectra/
    src/main_splines.cpp

## Interfaces and Dependencies

At the end of the refactor, these interfaces should exist:

In `include/experiments/shift_vertical_horizontal.h`, keep:

    void run_final_shift_tests(int nx, int ny, int horizontal_steps, int vertical_steps);

In `include/experiments/common_simulation_helpers.h`, define declarations for shared helpers. The exact set should be finalized as functions are moved, but the initial candidates are:

    Eigen::Matrix<double, 3, 2>
    calculateShapeDerivatives(const Eigen::Vector2d &p1,
                              const Eigen::Vector2d &p2,
                              const Eigen::Vector2d &p3);

    bool hasConnectivityChanged(const std::vector<ElementTriangle2D> &old_elements,
                                const std::vector<ElementTriangle2D> &new_elements,
                                const std::vector<size_t> &old_active,
                                const std::vector<size_t> &new_active);

    std::tuple<double, Eigen::Matrix2d, int> perform_remeshing_loop_reduction(...);

    void writeSizesToFile(int Nx, int Ny);

The `perform_remeshing_loop_reduction(...)` declaration must use the full existing signature from `src/main.cpp`. Do not simplify it during this refactor.

In `include/experiments/zanzotto_examples.h`, declare:

    void memory(int caller_id, int nx, int ny, int restart_iteration);
    void example_1_conti_zanzotto(int caller_id, int nx, int ny);
    void example_2_conti_zanzotto_triangular();

In `include/experiments/data_analysis.h`, declare:

    void analyze_data_from_folder(int caller_id, int nx, int ny,
                                  int iter_start, int iter_end, int n_eig);

In `include/experiments/deformation_examples.h`, declare:

    void debug_deformation_tests();
    void debug_deformation_tests_triangular();
    void example_1_shifting(int caller_id, int nx, int ny);

In `include/experiments/dislocation_indentation.h`, declare the functions that are intentionally callable from outside that module:

    void single_dislo_LJ();
    void indentation();

Keep `sortByCoordinates(...)`, `callback(...)`, and `savePositionEtaStress(...)` file-local if nothing else needs them. File-local means they should live in an unnamed namespace inside the `.cpp` file so their names do not leak across the program.

In `include/experiments/acoustic_studies.h`, declare:

    void parametricAcousticStudy();
    void parametricAcousticStudy_v2();

In `include/experiments/stress_controlled_examples.h`, declare:

    void example_3_stress_controlled_final_clean(int caller_id, int nx, int ny);

Keep dependencies local to each `.cpp` file. For example, if only acoustic studies need `include/acoustic_tensor.h`, include it in `src/experiments/acoustic_studies.cpp`, not in `src/main.cpp`.

## Change Log

- 2026-05-23 / Codex: Created the initial ExecPlan on branch `feature/main_factor` after reading `PLANS.md`, inspecting `src/main.cpp`, and recording the intended incremental refactor strategy. The reason for this change is to give a future implementer a safe, self-contained path for splitting `main.cpp` without changing calculations.

