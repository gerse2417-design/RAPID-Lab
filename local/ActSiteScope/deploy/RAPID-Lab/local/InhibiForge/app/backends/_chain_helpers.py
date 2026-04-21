"""
공유 헬퍼 — contigs 파싱 + protocol(binder/partial/fixbb) 감지.

run_mpnn.py / run_af2.py 양쪽에서 import 해서, ColabDesign의
designability_test.py에 있던 chain 처리 로직을 동일하게 재현한다.

원본:
  /home/sooyeon/amr/RFdiffusion/ColabDesign_repo/colabdesign/rf/designability_test.py
  (get_info: line 13-31, parse + protocol 분기: line 67-126)
"""
from string import ascii_uppercase, ascii_lowercase

import numpy as np


alphabet_list = list(ascii_uppercase + ascii_lowercase)


def get_info(contig):
    """
    한 contig 토큰을 (fixed_pos_list, [fixed_chain, free_chain])로 분해.

    예: "A1-80"     → ([1]*80, [True, False])    # 고정 모티프
        "100-100"  → ([0]*100, [False, True])   # 자유 설계 영역
        "A1-80/100-100" 같이 슬래시로 묶인 멀티 토큰은 미리 split 후 호출.
    """
    F = []
    free_chain = False
    fixed_chain = False
    sub_contigs = [x.split("-") for x in contig.split("/")]
    for n, sub in enumerate(sub_contigs):
        if len(sub) == 1:
            a = b = sub[0]
        else:
            a, b = sub[0], sub[1]
        if a[0].isalpha():
            L = int(b) - int(a[1:]) + 1
            F += [1] * L
            fixed_chain = True
        else:
            L = int(b)
            F += [0] * L
            free_chain = True
    return F, [fixed_chain, free_chain]


def parse_contigs(contigs_str):
    """
    '--contigs' 인자 문자열을 chain별 contig 리스트로 분해.

    구분자: ':' / ',' / 공백  (designability_test.py와 동일 규칙)
    각 contig 안의 '/0'은 ColabDesign 패딩 표현이라 제거.

    예: "A1-80/0 100-100" → ["A1-80", "100-100"]
    """
    contigs = []
    for contig_str in contigs_str.replace(" ", ":").replace(",", ":").split(":"):
        if len(contig_str) > 0:
            contig = []
            for x in contig_str.split("/"):
                if x != "0":
                    contig.append(x)
            contigs.append("/".join(contig))
    return contigs


def detect_protocol(contigs_str, *,
                    rm_aa="C",
                    copies=1,
                    use_multimer=False,
                    initial_guess=False,
                    param_dir="."):
    """
    contigs로부터 ColabDesign 프로토콜 자동 감지.

    Returns
    -------
    protocol : str  ("binder" | "partial" | "fixbb")
    af_model_kwargs : dict  →  mk_af_model(**af_model_kwargs)
    prep_flags : dict       →  af_model.prep_inputs(pdb_filename, **prep_flags)
    fixed_pos : list[int]   →  partial 프로토콜에서만 사용
                              (af_model.opt["fix_pos"] = np.where(...) 용)
    """
    contigs = parse_contigs(contigs_str)
    chains = alphabet_list[:len(contigs)]
    info = [get_info(x) for x in contigs]

    fixed_pos = []
    fixed_chains = []
    free_chains = []
    both_chains = []
    for pos, (fixed_chain, free_chain) in info:
        fixed_pos += pos
        fixed_chains += [fixed_chain and not free_chain]
        free_chains += [free_chain and not fixed_chain]
        both_chains += [fixed_chain and free_chain]

    common_flags = {
        "initial_guess": initial_guess,
        "best_metric": "rmsd",
        "use_multimer": use_multimer,
        "data_dir": param_dir,
        "model_names": ["model_1_multimer_v3" if use_multimer else "model_1_ptm"],
    }

    if sum(both_chains) == 0 and sum(fixed_chains) > 0 and sum(free_chains) > 0:
        # 고정 체인 + 자유 체인이 동시에 → binder 설계
        protocol = "binder"
        target_chains = [chains[n] for n, x in enumerate(fixed_chains) if x]
        binder_chains = [chains[n] for n, x in enumerate(fixed_chains) if not x]
        af_model_kwargs = {"protocol": "binder", **common_flags}
        prep_flags = {
            "target_chain": ",".join(target_chains),
            "binder_chain": ",".join(binder_chains),
            "rm_aa": rm_aa,
        }
    elif sum(fixed_pos) > 0:
        # 일부 잔기만 고정 → partial (fixbb + use_templates)
        protocol = "partial"
        af_model_kwargs = {"protocol": "fixbb", "use_templates": True, **common_flags}
        rm_template = np.array(fixed_pos) == 0
        prep_flags = {
            "chain": ",".join(chains),
            "rm_template": rm_template,
            "rm_template_seq": rm_template,
            "copies": copies,
            "homooligomer": copies > 1,
            "rm_aa": rm_aa,
        }
    else:
        # 전부 자유 → fixbb (백본 고정 시퀀스 설계)
        protocol = "fixbb"
        af_model_kwargs = {"protocol": "fixbb", **common_flags}
        prep_flags = {
            "chain": ",".join(chains),
            "copies": copies,
            "homooligomer": copies > 1,
            "rm_aa": rm_aa,
        }

    return protocol, af_model_kwargs, prep_flags, fixed_pos
