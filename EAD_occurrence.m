% This function is called to assess the presence of EADs.

function EAD_index = EAD_occurrence(time,Vm,tin)

tfin=tin+0.500e3; 
tin_roi=find(time>tin); tin_idx=tin_roi(1)-1;
tfin_roi=find(time>tfin); tfin_idx=tfin_roi(1);

time_roi = time(tin_idx:tfin_idx)-time(tin_idx);
Vm_roi = Vm(tin_idx:tfin_idx);

dVm_roi = (Vm_roi(2:end)-Vm_roi(1:end-1))./(time_roi(2:end)-time_roi(1:end-1));
dVm_dVm = dVm_roi(1:end-1).*dVm_roi(2:end);

[max_dVm idx_dVm] = max(dVm_roi);
dVm_0 = find(dVm_dVm(idx_dVm:end)<0);
dVm_0_count = length(dVm_0);

if dVm_0_count > 1,
    EAD_index = 1;
else
    EAD_index = 0;
end