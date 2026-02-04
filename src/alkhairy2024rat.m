addpath matlab_wrappers
pyenv(ExecutionMode="OutOfProcess")

function db = mag2db(x)
    db = 20*log10(x);
end

function fig1_2024rat()
    fils = [Filter(Ap=0.05, bp=1, Bu=2), Filter(Ap=0.05, bp=1, Bu=2.5), Filter(Ap=0.05, bp=1, Bu=3)];
    freqs = py.numpy.array(logspace(log10(0.8), log10(1.2), 10000));
    responses = arrayfun(@(fil) double(fil.PyFilter.filter.tf(1j*freqs)), fils, UniformOutput=false);

    tiledlayout(1, 2, TileSpacing='tight')

    nexttile
    twodb = mag2db(abs(responses{1}));
    semilogx(freqs, twodb-max(twodb))
    hold on
    twopointfivedb = mag2db(abs(responses{2}));
    semilogx(freqs, twopointfivedb-max(twopointfivedb))
    hold on
    threedb = mag2db(abs(responses{3}));
    semilogx(freqs, threedb-max(threedb))
    ylabel('Magnitude (dB)')
    xlabel('\beta')

    nexttile
    twocyc = unwrap(angle(responses{1})) / (2 * pi);
    semilogx(freqs, twocyc-max(twocyc))
    hold on
    twopointfivecyc = unwrap(angle(responses{2})) / (2 * pi);
    semilogx(freqs, twopointfivecyc-max(twopointfivecyc))
    hold on
    threecyc = unwrap(angle(responses{3})) / (2 * pi);
    semilogx(freqs, threecyc-max(threecyc))
    ylabel('Phase (cyc)')
    xlabel('\beta')
end

function fig2_2024rat()
    fils = arrayfun(@(Buval) Filter(Ap=0.05, bp=1, Bu=Buval), 2:8);
    % f = fils(1);
    % struct(f.get_computed_chars())
    Ncycs = arrayfun(@(fil) struct(fil.get_computed_chars()).Nbeta, fils);
    Qerbs = arrayfun(@(fil) struct(fil.get_computed_chars()).Qerb, fils);
    Q10s = arrayfun(@(fil) struct(fil.get_computed_chars()).Qn, fils);
    Q3s = arrayfun(@(fil) struct(fil.get_computed_chars()).Qn2, fils);
    Q15s = arrayfun(@(fil) struct(py.helpers.computedfiltercharacteristics(fil.PyFilter.filter.tf, n=15)).Qn, fils);

    plot(2:8, Ncycs, DisplayName='N_{cyc}')
    hold on
    plot(2:8, Qerbs, DisplayName='Q_{erb}')
    hold on
    plot(2:8, Q10s, DisplayName='Q_{10}')
    hold on
    plot(2:8, Q3s, DisplayName='Q_3')
    hold on
    plot(2:8, Q15s, DisplayName='Q_{15}')
    yline(0, ':k')
    legend
end

function fig3_2024rat()
    fil = Filter(Ap=0.05, bp=1, Bu=2, cf=1);
    function res = tones(t)
        res = 0;
        fis = [1, 5, 7/8, 1/5];
        tis = [20, 50, 70, 40];
        for idx = 1:4
            res = res + (exp(-((t-tis(idx))/5).^2) .* sin(2*pi*fis(idx)*t));
        end
    end
    fs = 100;
    timestamps = (1:(fs*100))/fs;
    sig = Signal(mode='t', data=tones(timestamps), fs=fs);
    outsig = fil.solve(sig, method='tf');

    tiledlayout(2, 1, TileSpacing='tight')

    nexttile
    plot(timestamps, sig.get_data('t'))
    ylabel('input')

    nexttile
    plot(timestamps, outsig.get_data('t'))
    ylabel('output')
end

function fig4_2024rat()
    fils1 = [Filter(Ap=0.15, bp=1, Bu=3), Filter(Ap=0.15, bp=1, Bu=5), Filter(Ap=0.045, bp=1, Bu=5)];
    fils2 = [Filter(Ap=0.15, bp=1, Bu=1.5), Filter(Ap=0.15, bp=1, Bu=2), Filter(Ap=0.15, bp=1, Bu=2.5), Filter(Ap=0.15, bp=1, Bu=3)];

    tiledlayout(2, 1, TileSpacing='tight')

    nexttile
    timestamps = 0:1999/10;
    for fil = fils1
        [~, ir] = fil.impulse_response_plot(times=timestamps, show=false);
        plot(timestamps, ir/max(ir))
        hold on
    end
    title('dependence of behavior of h on values of constants for integer B_u')
    xlabel('h($\tilde{t}$)', Interpreter='latex')
    ylabel('$\tilde{t}$', Interpreter='latex')
    yline(0, ':k')

    nexttile
    timestamps = 0:999/10;
    for fil = fils2
        [~, ir] = fil.impulse_response_plot(times=timestamps, show=false);
        plot(timestamps, ir/max(ir))
        hold on
    end
    title('h and phase of oscillatory component for integer and half-integer B_u')
    xlabel('h($\tilde{t}$)', Interpreter='latex')
    ylabel('$\tilde{t}$', Interpreter='latex')
    yline(0, ':k')
end

function fig7_2024rat()
    fil = Filter(Ap=0.1, bp=1, Bu=1.75, cf=1);
    sig = Signal(f_init=-2, f_final=2, fs=30, num_samples=3000);
    sig = Signal(mode='ttilde', data=sig.get_data('t'), fs=30);
    sol = fil.solve(sig, method='tf');

    tiledlayout(2, 2, TileSpacing='tight')

    nexttile
    plot(sig.PySignal.timestamps, sig.get_data('ttilde'))

    nexttile
    plot(sol.PySignal.timestamps, sol.get_data('ttilde'))

    nexttile
    % broken
end

% def fig7_2024rational():
%   fil = Filter(Ap=0.1, bp=1, Bu=1.75, cf=1)
%   sig = Signal.linear_chirp(f_init=-2, f_final=2, fs=30, num_samples=3000) # set this to be able to be ttilde
%   sig = Signal(mode='ttilde', data=sig.get_data('t'), fs=30)
%   sol = fil.solve(sig, method='tf')
% 
%   fig, axs = plt.subplots(2, 2, constrained_layout=True)
%   axs[0][0].plot(sig.timestamps, sig['ttilde'])
%   axs[0][1].plot(sol.timestamps, sol['ttilde'])
%   winsig, sigbound = sig.spectrogram(show=False)
%   winsol, solbound = sol.spectrogram(show=False)
%   axs[1][0].imshow(abs(winsig[0:21][::-1]), cmap='gray', aspect='auto', extent=(0, sigbound[1], 0, sigbound[3]/5))
%   axs[1][1].imshow(abs(winsol[0:21][::-1]), cmap='gray', aspect='auto', extent=(0, solbound[1], 0, solbound[3]/5))
%   plt.show()

function fig8_2024rat()
    function res = tones(t)
        res = 0;
        fis = [1, 5, 7/8, 1/5];
        tis = [20, 50, 70, 40];
        for idx = 1:4
            res = res + (exp(-((t-tis(idx))/5).^2) .* sin(2*pi*fis(idx)*t));
        end
    end
    fs = 10;
    ts = (0:(100*fs-1))/fs;
    data = tones(ts);
    sig = Signal(mode='ttilde', data=data, fs=fs);

    tiledlayout(2, 1, TileSpacing='tight')

    nexttile
    plot(ts, data)
    title('input')
    ylabel('$\tilde{t}$', Interpreter='latex')

    nexttile
    fil = Filter(Ap=0.15, bp=1, Bu=int32(5));
    for method = ["tf", "ir", "ode", "fde"]
        sol = fil.solve(sig, method=method).get_data('ttilde');
        plot(ts, sol/max(sol), ':', DisplayName=method);
        hold on
    end
    title('output')
    ylabel('$\tilde{t}$', Interpreter='latex')
    legend

    sgtitle('equivalence of representations A_p=0.15, b_p=1, B_u=5')
end

function fig9_2024rat()
    Ap = 0.04;
    bp = 1;
    Bu = 1.5;
    function res = myfunc(t)
        res = exp(-0.04.*t).*besselj(0, t);
    end
    ts = (0:999)/10;
    % Handles to MATLAB functions not handled?
    sig = Signal(mode='ttilde', data=myfunc(ts), fs=10);

    tiledlayout(2, 1, TileSpacing='tight')

    nexttile
    plot(sig.PySignal.timestamps, sig.get_data('ttilde'))
    title('input')

    nexttile
    fil = Filter(Ap=0.04, bp=1, Bu=1.5);
    for method = ["tf", "ir", "fde"]
        sol = fil.solve(sig, method=method).get_data('ttilde');
        plot(ts, sol/max(sol), ':', DisplayName=method);
        hold on
    end
    title('output')
    ylabel('$\tilde{t}$', Interpreter='latex')
    legend

    sgtitle('equivalence of representations A_p=0.04, b_p=1, B_u=1.5')
end

fig1_2024rat()
figure
fig2_2024rat()
figure
fig3_2024rat()
figure
fig4_2024rat()
% fig7_2024rat() % broken
figure
fig8_2024rat()
figure
fig9_2024rat()

% fil.PyFilter.filter.tf applied to an array errors. Redefine inside MATLAB
% objects? Force user to pass in numpy arrays

% char array vs string?