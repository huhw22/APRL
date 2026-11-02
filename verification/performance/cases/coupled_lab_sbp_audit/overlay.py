"""Materialize only read-only observer hooks on the frozen r5 source."""
from pathlib import Path
import shutil
def once(text,old,new):
    if text.count(old)!=1:raise RuntimeError('Observer overlay anchor mismatch: '+old[:100])
    return text.replace(old,new,1)
def apply(base:Path,target:Path,template:Path):
    shutil.copytree(base,target)
    for p in template.glob('*.h'):shutil.copy2(p,target/'src'/p.name)
    p=target/'src/lab_field_audit.h';s=p.read_text();s=once(s,'#include "paired_momentum_audit.h"','#include "paired_momentum_audit.h"\n#include "whole_cowan_budget.h"\n#include <cstdlib>')
    s=once(s,'std::unique_ptr<PairedMomentumAudit> paired_;','std::unique_ptr<PairedMomentumAudit> paired_;std::unique_ptr<WholeCowanBudget> whole_;')
    s=once(s,'const BoostFrameTransform& frame):g_','const BoostFrameTransform& frame,MPI_Comm comm=MPI_COMM_WORLD):g_')
    s=once(s,'    paired_.reset(new PairedMomentumAudit','    if(const char* path=std::getenv("APRL_WHOLE_COWAN_AUDIT")){if(!*path)throw std::runtime_error("Empty whole Cowan audit directory");whole_.reset(new WholeCowanBudget(path,rank,g,offset,globalNz,gx,gy,comm));}\n    paired_.reset(new PairedMomentumAudit')
    s=once(s,'oldB_[0]=f.bx();oldB_[1]=f.by();oldB_[2]=f.bz();sample','oldB_[0]=f.bx();oldB_[1]=f.by();oldB_[2]=f.bz();if(whole_)whole_->begin(f,step);sample')
    s=once(s,'  void magnetic(const EBFieldGrid& f,const EBElectricHaloView& h) {','  void magnetic(const EBFieldGrid& f,const EBElectricHaloView& h) {\n    if(whole_)whole_->magneticRaw(f,h);')
    s=once(s,'    paired_->sample(stage,f);','    if(whole_)whole_->sample(stage,f);\n    paired_->sample(stage,f);')
    s=once(s,'  void close(bool complete){paired_->close(complete);','  void beforeReferenceCurrent(const EBFieldGrid& f){if(whole_)whole_->beforeReferenceCurrent(f);}\n  void afterReferenceCurrent(const EBFieldGrid& f){if(whole_)whole_->afterReferenceCurrent(f);}\n  void close(bool complete){if(whole_)whole_->close(complete);paired_->close(complete);')
    p.write_text(s)
    p=target/'src/simulation.cpp';s=p.read_text()
    s=once(s,'    const char* labAuditDirectory=std::getenv("APRL_LAB_FIELD_AUDIT");','    const char* labAuditDirectory=std::getenv("APRL_LAB_FIELD_AUDIT");\n    if(std::getenv("APRL_WHOLE_COWAN_AUDIT")&&(!labAuditDirectory||!sourceCoupling_))throw std::invalid_argument("Whole Cowan audit requires the existing conservative-source lab audit");')
    s=once(s,'config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1],frame_));','config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1],frame_,communicator_));')
    old='        if(config_.mesh.movingWindow.referenceSurfaceCurrent)\n          applyMovingTransverseReferenceCurrent();';new='        if(labAudit)labAudit->beforeReferenceCurrent(*fields_);\n'+old+'\n        if(labAudit)labAudit->afterReferenceCurrent(*fields_);';s=once(s,old,new);p.write_text(s)
