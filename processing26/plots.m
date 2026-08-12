datadir = "E:\MATLAB-DATA2\kingair_data\test26\work\20260626_arr.c10.nc"
blurf = ncread(datadir,'alpha'); alpha=blurf(:);
blurf = ncread(datadir,'beta'); beta=blurf(:);
blurf = ncread(datadir,'tas'); tas=blurf(:);
blurf = ncread(datadir,'avroll'); avroll=blurf(:);
blurf = ncread(datadir,'avthead'); avthead=blurf(:);



ax1=subplot(5,1,1)
plot(alpha)
grid
v=axis;
axis([v(1) v(2) 0 15])
ylabel('alpha')

ax2=subplot(5,1,2)
plot(beta)
grid
v=axis;
axis([v(1) v(2) -5 5])
ylabel('beta')

ax3=subplot(5,1,3)
plot(tas)
grid
ylabel('tas')

ax4=subplot(5,1,4)
plot(avroll)
grid
ylabel('roll')

ax5=subplot(5,1,5)
plot(avthead)
grid
ylabel('thead')


axs = [ax1 ax2 ax3 ax4 ax5];

% Link x-limits
linkaxes(axs,'x');