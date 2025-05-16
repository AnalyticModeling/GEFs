addpath matlab_wrappers

initialRun

function res = filter_init()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    figure
    f1.bode_plot();

    % ir1 = pyrun('pyir = (lambda t: t*numpy.exp(-t)*numpy.sin(t))', 'pyir'); % -> currently invalid in MATLAB
    % f2 = Filter(ir=ir1);
    % figure
    % f2.bode_plot();

    f3 = Filter(coeffs={{1}, {1, 1, 1}});
    figure
    f3.bode_plot();

    f4 = Filter(roots={{1}, {1+2j, 1-2j}});
    figure
    f4.bode_plot();

    f5 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    figure
    f5.bode_plot();

    f6 = Filter(Bpeak=1.5, Nf=1.11, phiaccum=3.5, cf=10);
    figure
    f6.bode_plot();

    f7 = Filter(Ap=0.1, bp=1, Bu=3);
    figure
    f7.bode_plot();
end

function res = filter_multiband_params()
    f = Filter(Ap=0.1, bp=[0.5, 1, 1.5], Bu=[3, 5, 7], peak_magndb=1);
    f.bode_plot();
end

function res = filter_multiband_chars()
    f = Filter(Bpeak=[0.5, 1, 1.5], Nbeta=[15, 10, 5], phiaccum=3.5);
    f.bode_plot();
end

function res = filter_get_computed_chars()
    f = Filter(Bpeak=1, Nbeta=11.1, phiaccum=3.5);
    res = f.get_computed_chars();
end

function res = filter_get_computed_unnormalized_chars()
    f = Filter(Bpeak=1, Nf=11.1, phiaccum=3.5, cf=1);
    res = f.get_computed_unnormalized_chars();
end

function res = filter_get_orig_chars()
    f = Filter(Bpeak=1, Nbeta=11.1, phiaccum=3.5);
    res = f.get_orig_chars();
end

function res = filter_get_params()
    f = Filter(Ap=0.1, bp=1, Bu=3);
    res = f.get_params();
end

function res = filter_solve()
    % why is it printing the mode
    f = Filter(Ap=0.1, bp=1.0, Bu=2, cf=1);
    sig = Signal(mode='ttilde', func=(@(t) exp(-1/15.*(t-20).^2) .* cos(t)), num_samples=2000, fs=20);
    figure
    sig.plot();
    figure
    
    anstf = f.solve(sig, method='tf');
    anstf = rdivide(anstf, max(anstf.get_data('ttilde')));
    % anstf.plot(custom_title='tf solve')
    % figure

    ansir = f.solve(sig, method='ir');
    ansir = rdivide(ansir, max(ansir.get_data('ttilde')));
    % ansir.plot(custom_title='ir solve')
    % figure

    ansode = f.solve(sig, method='ode');
    ansode = rdivide(ansode, max(ansode.get_data('ttilde')));
    % ansode.plot(custom_title='ode solve')
    % figure

    ansfde = f.solve(sig, method='fde');
    ansfde = rdivide(ansfde, max(ansfde.get_data('ttilde')));
    % ansfde.plot(custom_title='fde solve')
    % figure

    plot(anstf.PySignal.timestamps, anstf.get_data('ttilde'), DisplayName='tf');
    hold on
    plot(ansir.PySignal.timestamps, ansir.get_data('ttilde'), '--', DisplayName='ir')
    hold on
    plot(ansode.PySignal.timestamps, ansode.get_data('ttilde'), DisplayName='ode')
    hold on
    plot(ansfde.PySignal.timestamps, ansfde.get_data('ttilde'), '--', DisplayName='fde')
    legend

    xlabel('Time (ms)')
end

function res = filter_bode()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.bode_plot();
    figure

    f2 = Filter(coeffs={{1}, {1, 1/2, 1/4}});
    f2.bode_plot();
    figure

    f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f3.bode_plot();
end

function res = filter_frequency_real_imag()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.frequency_real_imag_plot();
    figure

    f2 = Filter(coeffs={{1}, {1, 1/2, 1/4}});
    f2.frequency_real_imag_plot();
    figure

    f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f3.frequency_real_imag_plot();
end

function res = filter_nichols()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.nichols_plot();
    figure

    f2 = Filter(coeffs={{1}, {1, 1/2, 1/4}});
    f2.nichols_plot();
    figure

    f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f3.nichols_plot();
end

function res = filter_nyquist()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.nyquist_plot();
    figure

    f2 = Filter(coeffs={{1}, {1, 1/2, 1/4}});
    f2.nyquist_plot();
    figure

    f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f3.nyquist_plot();
end

function res = filter_ir()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.impulse_response_plot(times=linspace(0, 10, 300));
    figure

    f2 = Filter(coeffs={{1}, {1, 1/2, 1/4}});
    f2.impulse_response_plot(times=linspace(0, 10, 300));
    figure

    f3 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f3.impulse_response_plot(times=linspace(0, 100, 1000));
    figure

    f4 = Filter(Ap=0.01, bp=1.0, Bu=3, cf=1/2/pi);
    f4.impulse_response_plot(times=linspace(0, 100, 1000));
end

function res = filter_pz()
    f1 = Filter(type='V', Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f1.pzplot();
    figure

    f2 = Filter(coeffs={{1, 2}, {1, 1/2, 1/4}});
    f2.pzplot();
end

function res = filter_Qns()
    fil = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    fil.Qn_plot();
end

function res = filter_characteristic_error()
    fil = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    fil.characteristic_error();
end

function res = signal_init()
    ts1 = 0:99;
    s1 = Signal(mode='t', data=(1+ts1/100).*sin(ts1/10));
    s1.plot()
    figure

    s2 = Signal(mode='f', data=(0:19)/10);
    s2.plot()
    figure

    s3 = Signal(mode='t', func=(@(x) (1+x/100).*sin(x/10)), num_samples=100);
    s3.plot()
    figure

    s4 = Signal(mode='f', func=(@(x) x/10), num_samples=20);
    s4.plot()
    figure

    s5 = Signal(f_init=1.5, f_final=6, fs=100, num_samples=800);
    s5.plot()
    figure

    s6 = Signal(freq_func=(@(x) x/10), fs=20, num_samples=200);
    s6.plot()
end

function res = signal_arith()
    sig = Signal(f_init=1.5, f_final=6, fs=100, num_samples=800);
    sig.plot()
    figure

    newsig1 = sig+1;
    newsig1.plot()
    figure

    newsig2 = sig.*sig;
    newsig2.plot()
end

function res = signal_at_time()
    sig = Signal(func=(@(x) sin(x/1.3)), fs=10, num_samples=100);
    sig.plot()
    sig.at_time(8.1) % -0.0524
    sig.at_time(2.6*pi) % -0.0068 (2.6*pi is about 8.16814)
    sig.at_time(8.2) % 0.0245
end

function [res1, res2] = signal_get_data()
    sig = Signal(f_init=1.5, f_final=6, fs=100, num_samples=800);
    res1 = sig.get_data('t');
    res2 = sig.get_data('f');
end

function res = signal_resample()
    sig = Signal(f_init=1.5, f_final=6, fs=100, num_samples=800);
    sig.plot()
    figure
    
    newsig = sig.resample(77);
    newsig.plot()
end

function res = signal_envelope_analytic()
    sig = Signal(f_init=1.5, f_final=6, fs=100, num_samples=800);
    newsig = 1 + sig.*(1.5+sin((0:799)/65));
    xaxis = double(0:newsig.length-1)/double(100.0);
    [upper, lower] = newsig.envelope_analytic();
    plot(xaxis, newsig.get_data('t'));
    hold on
    plot(xaxis, upper);
    hold on
    plot(xaxis, lower);
    xlabel('Time (ms)')
end

function res = filterbank_add()
    % Conversion of MATLAB 'Filter' to Python is not supported.
    % does this mean rewrite

    % fs = FilterBank(topology="parallel", Ap=[0.1, 0.1], bp=[0.5, 1.0], Bu=[3, 3]);
    % fil = Filter(Ap=0.1, bp=1.5, Bu=3).PyFilter
    % fs.add(fil, source=fs.filters(end));
    % Filter(Ap=0.1, bp=1.5, Bu=3).PyFilter
    % fs.bode_plot()
end

function res = filterbank_process_signal()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    os = fs.process_signal(Signal(f_init=1, f_final=10, fs=100, num_samples=200));

    for sig = os.outSignals
        figure
        sig.plot()
    end
end

function res = filterbank_bode()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    fs.bode_plot()
end

function res = outputsignals_init_kindof()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    insig = Signal(f_init=1, f_final=10, fs=100, num_samples=200);
    insig.plot()
    os = fs.process_signal(insig);

    for sig = os.outSignals
        figure
        sig.plot()
    end
end

function res = outputsignals_autocorrelates()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    os = fs.process_signal(Signal(f_init=1, f_final=10, fs=100, num_samples=200));
    os.autocorrelates()
end

function res = outputsignals_correlate_with()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    os = fs.process_signal(Signal(f_init=1, f_final=10, fs=100, num_samples=200));
    sig = Signal(f_init=2, f_final=5, fs=100, num_samples=200);
    figure
    sig.plot()
    figure
    os.correlate_with(sig)
end

function res = outputsignals_correlogram()
    fs = FilterBank(topology="parallel", Ap=[0.1, 0.1, 0.1], bp=[0.5, 1.0, 1.5], Bu=[3, 3, 3]);
    os = fs.process_signal(Signal(f_init=1, f_final=10, fs=100, num_samples=200));
    os.correlogram()
end

function res = cochlea_init()
    c = Cochlea(Ap=0.3768*exp(-0.1*(0:3)), bp=[0.5, 1, 1.5, 2], Bu=3.714*exp(0.03*(0:3)));
    res = c.bode_plot();
end

function res = cochlea_at_location()
    c = Cochlea(Ap=0.3768*exp(-0.1*(0:3)), bp=[0.5, 1, 1.5, 2], Bu=3.714*exp(0.03*(0:3)));
    fil = c.filter_at_location(0);
    res = fil.bode_plot();
end

function res = cochlea_wavenumber()
    c = Cochlea(Ap=0.3768*exp(-0.1*(0:3)), bp=[0.5, 1, 1.5, 2], Bu=3.714*exp(0.03*(0:3)));
    res = c.plot_wavenumber();
end

function res = cochlea_impedance()
    c = Cochlea(Ap=0.3768*exp(-0.1*(0:3)), bp=[0.5, 1, 1.5, 2], Bu=3.714*exp(0.03*(0:3)));
    res = c.plot_impedance();
end

function res = cochlea_heatmap()
    c = Cochlea(type='V', aAp=0.3768, bAp=-0.1366, bp=[1, 1, 1, 1, 1], aBu=3.714, bBu=0.03123, xs=[0, 1, 2, 3, 4]);
    function innerres = tones(t)
        innerres = 0;
        fis = [1.5, 8, 1.5, 3];
        tis = [200, 400, 700, 400];
        for idx = 1:4
            innerres = innerres + (exp(-((t-tis(idx))/50).^2) .* sin(2*pi*fis(idx)*t));
        end
    end
    sample_rate = 1;
    % endtime = 1000*sample_rate;
    endtime = 100;
    timestamps = 0:(endtime-1);
    sig = Signal(mode='t', data=timestamps./sample_rate, fs=sample_rate);
    res = c.signal_response_heatmap(sig);
    % res{0}
end

% 
% def cochlea_heatmap():
%   c = Cochlea.five_param(type='V', aAp=0.3768, bAp=-0.1366, bp=[1, 1, 1, 1, 1], aBu=3.714, bBu=0.03123, xs=[i for i in range(5)])
%   pairs = [(1.5, 200), (8, 400), (1.5, 700), (0.3, 400)]
%   def tones(t):
%     ans = 0
%     for i in range(4):
%       fi, ti = pairs[i]
%       ans += np.exp(-((t-ti)/50)**2) * np.sin(2*np.pi*fi*t)
%     return ans
%   sig = Signal(mode='t', data=[tones(t/100) for t in range(100000)], fs=100)
%   c.signal_response_heatmap(sig)


% filter_init()
% filter_multiband_params()
% filter_multiband_chars()
% filter_get_computed_chars()
% filter_get_computed_unnormalized_chars()
% filter_get_orig_chars()
% filter_get_params()
% filter_solve()
% filter_bode()
% filter_frequency_real_imag()
% filter_nichols()
% filter_nyquist()
% filter_ir()
% filter_pz()
% filter_Qns()
% filter_characteristic_error()
% signal_init()
% signal_arith()
% signal_at_time()
% signal_get_data()
% signal_resample()
% signal_envelope_analytic()
% filterbank_add()
% filterbank_process_signal()
% filterbank_bode()
% outputsignals_init_kindof()
% outputsignals_autocorrelates()
% outputsignal_correlate_with()
outputsignals_correlogram()
% cochlea_init()
% cochlea_at_location()
% cochlea_wavenumber()
% cochlea_impedance()
% cochlea_heatmap()

function res = filter_test()
    % gcf;
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    figure
    f1.bode_plot();
end

function res = test2()
    f = Filter(Ap=0.1, bp=[0.5, 1, 1.5], Bu=[3, 5, 7], peak_magndb=1);
    [freqs, magns, phases] = f.bode_plot();
    plot(freqs, magns)
end