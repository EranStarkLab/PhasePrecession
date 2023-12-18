
load('1_percession_only.mat');
spk_phs_spec(spk, phs, periods );
sgtitle('Unit 1: precession only');

load('2_lock_only.mat');
spk_phs_spec(spk, phs, periods );
sgtitle('Unit 2: lock only');

load('3_percession_lock.mat')
spk_phs_spec(spk, phs, periods );
sgtitle('Unit 3: lock and precession');