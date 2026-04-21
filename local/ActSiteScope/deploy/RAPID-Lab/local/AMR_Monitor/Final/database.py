import os
from sqlalchemy import create_engine, Column, Integer, String, Float, JSON, ForeignKey, DateTime, Boolean
from sqlalchemy.orm import declarative_base, sessionmaker, relationship
from datetime import datetime

Base = declarative_base()

class PipelineRunDB(Base):
    """
    파이널 워크플로우를 관통하는 메인 레코드 테이블
    전체 실행 정보(단백질 기본 정보 및 파이프라인 모드)를 담습니다.
    """
    __tablename__ = 'pipeline_run_db'

    run_id = Column(String, primary_key=True)         # ex: TARGET_001_run_1
    mode = Column(String, nullable=False)             # "RADAR_BOARDS" or "Full_Pipeline"
    target_fasta_path = Column(String, nullable=True) # Pipeline 모드일때 단백질 시퀀스
    created_at = Column(DateTime, default=datetime.utcnow)

    # 역관계: Stage 1 ~ 최종 결과
    stage1_targets = relationship("TargetStructureDB", back_populates="run", cascade="all, delete-orphan")
    radar_boards = relationship("RadarBoardsDB", back_populates="run", cascade="all, delete-orphan", uselist=False)


class RadarBoardsDB(Base):
    """
    RADAR -> BOARDS DB 모드 선택 시 사용되는 데이터베이스
    """
    __tablename__ = 'radar_boards_db'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(String, ForeignKey('pipeline_run_db.run_id'), nullable=False)
    radar_output = Column(JSON, nullable=True) # RADAR 분석 결과 (내성/기능 예측)
    boards_entry_id = Column(String, nullable=True) # BOARDS DB 매칭된 ID
    
    run = relationship("PipelineRunDB", back_populates="radar_boards")


class TargetStructureDB(Base):
    """
    [Stage 1] RADAR / AlphaFold2 - AMR 유전자 탐지 및 항생제 내성 유발 단백질 3D 구조 예측
    """
    __tablename__ = 'target_structure_db'

    model_id = Column(String, primary_key=True) # run_id + "_rank_1" 형태
    run_id = Column(String, ForeignKey('pipeline_run_db.run_id'), nullable=False)
    
    pdb_path = Column(String, nullable=True)
    plddt_score = Column(Float, nullable=True)
    pae_data_path = Column(String, nullable=True)
    
    run = relationship("PipelineRunDB", back_populates="stage1_targets")
    hotspots = relationship("HotspotCandidateDB", back_populates="target_model", cascade="all, delete-orphan")


class HotspotCandidateDB(Base):
    """
    [Stage 2] 단백질 작용 부위 예측
    """
    __tablename__ = 'hotspot_candidate_db'

    id = Column(Integer, primary_key=True, autoincrement=True)
    model_id = Column(String, ForeignKey('target_structure_db.model_id'), nullable=False)
    
    hotspot_residues = Column(JSON, nullable=True)
    score_data = Column(JSON, nullable=True)
    rfdiffusion_input_path = Column(String, nullable=True) 

    target_model = relationship("TargetStructureDB", back_populates="hotspots")
    designed_binders = relationship("BinderDesignDB", back_populates="hotspot", cascade="all, delete-orphan")
    

class BinderDesignDB(Base):
    """
    [Stage 3 & Stage 4] 향후 연동 준비용
    Stage 3: RFdiffusion - 억제 단백질 3D 구조 설계 (De novo)
    Stage 4: ProteinMPNN - 억제 단백질 아미노산 서열 생성
    """
    __tablename__ = 'binder_design_db'
    
    binder_id = Column(String, primary_key=True) # hotspot_id + "_binder_1"
    hotspot_id = Column(Integer, ForeignKey('hotspot_candidate_db.id'), nullable=False)
    
    rfdiffusion_pdb_path = Column(String, nullable=True)   # 생성된 백본 구조
    mpnn_fasta_path = Column(String, nullable=True)        # 생성된 바인더 아미노산 서열 (여러개일 수 있으나 여기선 1:1 매핑)
    mpnn_score = Column(Float, nullable=True)              # 설계 적합성 점수
    
    hotspot = relationship("HotspotCandidateDB", back_populates="designed_binders")
    docking_evals = relationship("DockingEvaluationDB", back_populates="binder", cascade="all, delete-orphan")
    dna_works = relationship("DNAWorksDB", back_populates="binder", cascade="all, delete-orphan", uselist=False)


class DockingEvaluationDB(Base):
    """
    [Stage 5] 향후 연동 준비용
    AlphaFold-Multimer + LightDock - 결합체 구조 예측 및 검증
    """
    __tablename__ = 'docking_evaluation_db'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    binder_id = Column(String, ForeignKey('binder_design_db.binder_id'), nullable=False)
    
    complex_pdb_path = Column(String, nullable=True)       # 결합 상태의 구조
    af_multimer_plddt = Column(Float, nullable=True)       # pLDDT or ipTM 
    lightdock_score = Column(Float, nullable=True)
    dockq_score = Column(Float, nullable=True)             # DockQ 점수
    
    binder = relationship("BinderDesignDB", back_populates="docking_evals")


class DNAWorksDB(Base):
    """
    [Stage 6] 향후 연동 준비용
    DNAWorks - DNA 합성 설계 (올리고뉴클레오타이드)
    """
    __tablename__ = 'dnaworks_db'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    binder_id = Column(String, ForeignKey('binder_design_db.binder_id'), nullable=False)
    
    dna_sequence = Column(String, nullable=True)
    optimization_score = Column(Float, nullable=True)
    
    binder = relationship("BinderDesignDB", back_populates="dna_works")


def init_db(db_path="sqlite:///final_workflow.db"):
    """Initialize the SQLite database for the Final Workflow."""
    engine = create_engine(db_path, echo=False)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()

if __name__ == "__main__":
    session = init_db()
    print("Database initialized successfully at final_workflow.db")
    session.close()
