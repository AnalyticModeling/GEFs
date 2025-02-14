addpath matlab_wrappers
pyenv(ExecutionMode="OutOfProcess")

function fig2_2022
    tiledlayout(2, 2, TileSpacing='tight')
    
    % first two cochleas
    c1 = Cochlea(Ap=[0.11], bp=[1], Bu=[7], CF0=1);
    [betas1, reals1, imags1] = c1.plot_wavenumber(betas=logspace(log10(0.5), log10(2), 10000), show=false);
    c2 = Cochlea(Ap=[0.055], bp=[1], Bu=[7], CF0=10);
    [betas2, reals2, imags2] = c2.plot_wavenumber(betas=logspace(log10(0.5), log10(2), 10000), show=false);

    nexttile
    semilogx(betas1, reals1)
    hold on
    semilogx(betas2, reals2, '--')
    ylabel('Re\{k\} (mm^{-1})')
    xlabel('\beta')
    yline(0, ':k')
    xline(1, ':k')
    
    nexttile
    semilogx(betas1, imags1)
    hold on
    semilogx(betas2, imags2, '--')
    ylabel('Im\{k\} (mm^{-1})')
    xlabel('\beta')
    yline(0, ':k')
    xline(1, ':k')
    
    % last two cochleas
    c1 = Cochlea(Ap=[0.11], bp=[1], Bu=[7], CF0=1);
    [betas1, reals1, imags1] = c1.plot_impedance(betas=logspace(log10(0.5), log10(2), 10000), show=false);
    c2 = Cochlea(Ap=[0.055], bp=[1], Bu=[7], CF0=10);
    [betas2, reals2, imags2] = c2.plot_impedance(betas=logspace(log10(0.5), log10(2), 10000), show=false);

    nexttile
    semilogx(betas1, reals1)
    hold on
    semilogx(betas2, reals2, '--')
    ylabel('Re\{Z\}/(2\pi\rho l)')
    xlabel('\beta')
    yline(0, ':k')
    xline(1, ':k')
    
    nexttile
    semilogx(betas1, imags1)
    hold on
    semilogx(betas2, imags2, '--')
    ylabel('Im\{Z\}/(2\pi\rho l)')
    xlabel('\beta')
    yline(0, ':k')
    xline(1, ':k')
end

fig2_2022