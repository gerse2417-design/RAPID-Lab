"""
amr_dock_backend/__init__.py
──────────────────────────────────────────────────────────────
AMR 도킹 백엔드 패키지 진입점.

이 파일은 모든 서브모듈의 public 심볼을 단일 네임스페이스로 재노출한다.
프론트엔드(app_amr_dock_interactive_v12.py)는 이 파일 하나에서 모든 것을 import한다:

    from amr_dock_backend import (
        parse_rank_list, run_full_docking_pipeline, get_tool_status, ...
    )

책임 범위:
  - 서브모듈 간 의존성를 캡슐화하고 public API 경계를 명확히 함
  - 기존 run_amr_dock_interactive_v10.py 에서 import 하던 심볼 모두 재노출 (하위 호환)

의존성:
  - 내부: amr_dock_backend 의 모든 서브모듈
"""

from .pdb_utils import (
    get_sequence_from_pdb,
    preprocess_pdb_for_docking,
    clean_pdb_for_lightdock,
    substitute_nonstandard_residues,
    run_pdbfixer_pipeline,
    NONSTANDARD_RESIDUE_MAP,
    PDBFIXER_AVAILABLE,
)

from .parsers import (
    parse_rank_list,
    parse_cluster_repr,
    parse_rank_filtered,
    parse_dnaworks_logfile,
)

from .swarm_utils import (
    iter_swarm_dirs,
    collect_all_swarm_clusters,
    get_swarm_pdb_files,
)

from .restraints import (
    parse_atom_input,
    parse_pdb_atoms,
    generate_restraints_list,
)

from .scripts import (
    build_pymol_heatmap_script,
    build_chimerax_script,
    build_dnaworks_input,
)

from .analysis import (
    find_interface_residues,
    build_energy_landscape_df,
    parse_swarm_coordinates,
)

from .env import (
    get_tool_status,
    get_optimal_cores,
    DNAWORKS_EXECUTABLE,
    CHIMERAX_BIN,
)

from .subprocess_utils import (
    run_subprocess,
    save_timing,
    format_duration_min,
)

from .pipeline import (
    run_full_docking_pipeline,
    run_bsas_clustering,
)

from .job import get_job_summary

__all__ = [
    # pdb_utils
    "get_sequence_from_pdb", "preprocess_pdb_for_docking",
    "clean_pdb_for_lightdock", "substitute_nonstandard_residues",
    "run_pdbfixer_pipeline", "NONSTANDARD_RESIDUE_MAP", "PDBFIXER_AVAILABLE",
    # parsers
    "parse_rank_list", "parse_cluster_repr", "parse_rank_filtered",
    "parse_dnaworks_logfile",
    # swarm_utils
    "iter_swarm_dirs", "collect_all_swarm_clusters", "get_swarm_pdb_files",
    # restraints
    "parse_atom_input", "parse_pdb_atoms", "generate_restraints_list",
    # scripts
    "build_pymol_heatmap_script", "build_chimerax_script", "build_dnaworks_input",
    # analysis
    "find_interface_residues", "build_energy_landscape_df", "parse_swarm_coordinates",
    # env
    "get_tool_status", "get_optimal_cores", "DNAWORKS_EXECUTABLE", "CHIMERAX_BIN",
    # subprocess_utils
    "run_subprocess", "save_timing", "format_duration_min",
    # pipeline
    "run_full_docking_pipeline", "run_bsas_clustering",
    # job
    "get_job_summary",
]
