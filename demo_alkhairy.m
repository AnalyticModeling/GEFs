addpath matlab_wrappers

function res = filter_init()
    tf1 = pyrun('pytf = (lambda s: 1/(1+s+s**2))', 'pytf');
    f1 = Filter(tf=tf1);
    f1.bode_plot();
    figure

    % ir1 = pyrun('pyir = (lambda t: np.sin(t)/t if t != 0 else 1)', 'pyir');
    % f2 = Filter(ir=ir1);
    % f2.bode_plot();
    % figure

    f3 = Filter(coeffs={{1}, {1, 1, 1}});
    f3.bode_plot();
    figure

    f4 = Filter(roots={{1}, {1+2j, 1-2j}});
    f4.bode_plot();
    figure

    f5 = Filter(Bpeak=1.5, Nbeta=11.1, phiaccum=3.5);
    f5.bode_plot();
    figure

    f6 = Filter(Bpeak=1.5, Nf=1.11, phiaccum=3.5, cf=10);
    f6.bode_plot();
    figure

    f7 = Filter(Ap=0.1, bp=1, Bu=3);
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
    
    anstf = f.solve(sig, method='tf');
    anstf = rdivide(anstf, max(anstf.get_data('ttilde')));
    anstf.plot(custom_title='tf solve')
    figure

    ansir = f.solve(sig, method='ir');
    ansir = rdivide(ansir, max(ansir.get_data('ttilde')));
    ansir.PySignal.mode
    ansir.plot(custom_title='ir solve')
    figure

    ansode = f.solve(sig, method='ode');
    ansode = rdivide(ansode, max(ansode.get_data('ttilde')));
    ansode.plot(custom_title='ode solve')
    figure

    ansfde = f.solve(sig, method='fde');
    ansfde = rdivide(ansfde, max(ansfde.get_data('ttilde')));
    ansfde.plot(custom_title='fde solve')
    figure

    plot(anstf.PySignal.timestamps, anstf.get_data('ttilde'), DisplayName='tf');
    hold on
    plot(ansir.PySignal.timestamps, ansir.get_data('ttilde'), '--', DisplayName='ir')
    hold on
    plot(ansode.PySignal.timestamps, ansode.get_data('ttilde'), DisplayName='ode')
    hold on
    plot(ansfde.PySignal.timestamps, ansfde.get_data('ttilde'), '--', DisplayName='fde')
    legend
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
filter_characteristic_error()