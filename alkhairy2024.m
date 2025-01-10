addpath matlab_wrappers
pyenv(ExecutionMode="OutOfProcess")

function db = mag2db(x)
    db = 20*log10(x);
end

function fig1_2024
    Ap = 0.055;
    bp = 1;
    Bu = 7;
    p = 1j*bp - Ap;
    P = @(s) ((s - p) .* (s - conj(p))).^(-Bu);
    Psharp = @(s) (s-p).^(-Bu);
    k = @(beta) Bu * (1./(1j*beta - p) + 1./(1j*beta - conj(p)));
    ksharp = @(beta) Bu * (1./(1j*beta - p));
    freqs = logspace(log10(0.1), log10(2), 10000);

    respP = P(1j*freqs);
    respPsharp = Psharp(1j*freqs);
    respk = k(freqs);
    respksharp = ksharp(freqs);

    tiledlayout(2, 2, TileSpacing='tight')

    nexttile
    Pdb = mag2db(abs(respP));
    Psharpdb = mag2db(abs(respPsharp));
    semilogx(freqs, Pdb-max(Pdb))
    hold on
    semilogx(freqs, Psharpdb-max(Psharpdb), '--')
    title('mag\{P\} (dB)')
    xlabel('\beta')
    
    nexttile
    kre = real(respk);
    ksharpre = real(respksharp);
    semilogx(freqs, kre)
    hold on
    semilogx(freqs, ksharpre, '--')
    title('Re\{k_\beta\}')
    xlabel('\beta')

    nexttile
    Pcyc = unwrap(angle(respP)) / (2 * pi);
    Psharpcyc = unwrap(angle(respPsharp)) / (2 * pi);
    semilogx(freqs, Pcyc-Pcyc(1))
    hold on
    semilogx(freqs, Psharpcyc-Psharpcyc(1), '--');
    title('phase\{P\} (cyc)')
    xlabel('\beta')
    
    nexttile
    kim = imag(respk);
    ksharpim = imag(respksharp);
    semilogx(freqs, kim)
    hold on
    semilogx(freqs, ksharpim, '--')
    title('Re\{k_\beta\}')
    xlabel('\beta')
end

function fig3_2024
    filP = Filter(Bpeak=1, Nbeta=19.1, Qerb=25.9);
    filV = Filter(type='V', Bpeak=1, Nbeta=19.1, Qerb=25.9);
    
    P = filP.PyFilter.filter.tf;
    Psharpoutput = cell(filP.PyFilter.filter.sharp_approximation());
    Psharp = Psharpoutput{1};
    V = filV.PyFilter.filter.tf;
    freqs = py.numpy.array(logspace(log10(0.1), log10(8), 10000));

    respP = double(P(1j*freqs));
    respPsharp = double(Psharp(1j*freqs));
    respV = double(V(1j*freqs));

    tiledlayout(1, 2, TileSpacing='tight')

    nexttile
    Pdb = mag2db(abs(respP));
    Psharpdb = mag2db(abs(respPsharp));
    Vdb = mag2db(abs(respV));
    semilogx(freqs, Pdb-max(Pdb))
    hold on
    semilogx(freqs, Psharpdb-max(Psharpdb), '--')
    hold on
    semilogx(freqs, Vdb-max(Vdb), ':')
    ylabel('Magnitude (dB)')
    xlabel('\beta')
    
    nexttile
    Pcyc = unwrap(angle(respP)) / (2 * pi);
    Psharpcyc = unwrap(angle(respPsharp)) / (2 * pi);
    Vcyc = unwrap(angle(respV)) / (2 * pi);
    semilogx(freqs, Pcyc-Pcyc(1))
    hold on
    semilogx(freqs, Psharpcyc-Psharpcyc(1), '--')
    hold on
    semilogx(freqs, Vcyc-Vcyc(1), ':')
    ylabel('Phase (cyc)')
    xlabel('\beta')
end

function fig5_2024
    f1 = Filter(Bpeak=1, Qn=7.8, Sbeta=1.7e3, n=3);
    f2 = Filter(Bpeak=1.5, Qn=12, Sbeta=1.7e3, n=3);
    f3 = Filter(Bpeak=2, Qn=16, Sbeta=1.7e3, n=3);
    fall = Filter(pyfilter=py.Filter.Filter.multiband_chars(Bpeak=[1, 1.5, 2], Qn=[7.8, 12, 16], Sbeta=1.7e3, n=3));
    betas = linspace(0.1, 2.5, 10000);
    
    minidx1 = 4791;
    minidx2 = 6874;
    [f1f, f1m, f1p] = f1.bode_plot(freqs=betas, show=false);
    [f2f, f2m, f2p] = f2.bode_plot(freqs=betas, show=false);
    [f3f, f3m, f3p] = f3.bode_plot(freqs=betas, show=false);
    [fallf, fallm, fallp] = fall.bode_plot(freqs=betas, show=false);
    % data = cellfun(@(fil) [fil.bode_plot(freqs=betas, show=false)], {f1, f2, f3, fall}, UniformOutput=false);
    tiledlayout(2, 2, TileSpacing='tight')

    nexttile
    plot(f1f, f1m)
    hold on
    plot(f2f, f2m)
    hold on
    plot(f3f, f3m)
    hold on
    plot(fallf, fallm)

    nexttile
    plot(f1f, f1p)
    hold on
    plot(f2f, f2p)
    hold on
    plot(f3f, f3p)
    hold on
    plot(fallf, fallp)

    nexttile([1, 2])
    c1 = struct(py.helpers.computedfiltercharacteristics(tfunc=fall.PyFilter.filter.tf, betas=betas(1:minidx1), n=3));
    c2 = struct(py.helpers.computedfiltercharacteristics(tfunc=fall.PyFilter.filter.tf, betas=betas(minidx1+1:minidx2), n=3));
    c3 = struct(py.helpers.computedfiltercharacteristics(tfunc=fall.PyFilter.filter.tf, betas=betas(minidx2+1:end), n=3));
    
    labelcell = {'1:\beta_{peak}', '2:\beta_{peak}', '3:\beta_{peak}', '1:Q_3', '2:Q_3', '3:Q_3', '1:S_\beta', '2:S_\beta', '3:S_\beta'};
    labels = categorical(labelcell);
    labels = reordercats(labels, labelcell);
    expected = [1, 1.5, 2, 7.8, 12, 16, 1.7e3, 1.7e3, 1.7e3];
    estimated = [c1.Bpeak, c2.Bpeak, c3.Bpeak, c1.Qn, c2.Qn, c3.Qn, c1.Sbeta, c2.Sbeta, c3.Sbeta];
    errors = zeros(1, 9);
    for i = 1:9
        errors(i) = abs(estimated(i)/expected(i)-1);
    end
    bar(labels, errors)
end

% fig1_2024
% fig3_2024
% fig5_2024

% fil.PyFilter.filter.tf applied to an array errors. Redefine inside MATLAB
% objects? Force user to pass in numpy arrays