addpath matlab_wrappers
pyenv(ExecutionMode="OutOfProcess")
% axis labeling?

function fig3_2019
    c = Cochlea(Ap=[0.05], bp=[1], Bu=[1.3]);
    [betas, reals, imags] = c.plot_wavenumber(betas=logspace(log10(0.7), log10(1.2), 10000), show=false);
    tiledlayout(2, 1, TileSpacing='tight')
    
    nexttile
    semilogx(betas, reals, '--')
    yline(0, ':k')
    xline(1, ':k')
    xlabel('\beta')
    ylabel('Re\{k\} (mm^{-1})')
    
    nexttile
    semilogx(betas, imags, '--')
    yline(0, ':k')
    xline(1, ':k')
    xlabel('\beta')
    ylabel('Im\{k\} (mm^{-1})')
end

function fig6_2019
    c = Cochlea(type='V', aAp=0.3768, bAp=-0.1366, bp=[0.2, 0.5, 2, 5, 15], aBu=3.714, bBu=0.03123, xs=[0, 1, 2, 3, 4]);
    fils = c.bode_plot(freqs=logspace(log10(0.1), log10(30), 10000), show=false);
    tiledlayout(2, 1, TileSpacing='tight')
    sgtitle('Q and N of velocities of a cochlea')
    
    nexttile
    for i = 1:5
        fil = fils{i};
        chars = c.filters(i).get_computed_chars();
        xaxis = fil{1};
        magn = fil{2};
        semilogx(xaxis, magn)
        hold on
        xlabel('Normalized frequency')
        ylabel('Magnitude (dB)')
    end
    
    nexttile
    for i = 1:5
        fil = fils{i};
        chars = c.filters(i).get_computed_chars();
        xaxis = fil{1};
        phase = fil{3};
        semilogx(xaxis, phase)
        hold on
        xlabel('Normalized frequency')
        ylabel('Phase (radians)')
    end
end

function fig8_2019
    function res = tones(t)
        res = 0;
        fis = [1.5, 8, 1.5, 3];
        tis = [200, 400, 700, 400];
        for idx = 1:4
            res = res + (exp(-((t-tis(idx))/50).^2) .* sin(2*pi*fis(idx)*t));
        end
    end

    c = Cochlea(type='V', aAp=0.3768, bAp=-0.1366, bp=[0.2, 0.5, 2, 5, 15], aBu=3.714, bBu=0.03123, xs=[0, 1, 2, 3, 4]);
    sig = Signal(mode='t', data=tones((1:10000)/10), fs=10);
    c.signal_response_heatmap(sig, len_xs=20)
end

fig3_2019
figure
fig6_2019
figure
fig8_2019