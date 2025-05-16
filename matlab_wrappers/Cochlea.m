classdef Cochlea < FilterBank
    properties
        PyCochlea
    end
    methods
        function obj = Cochlea(options)
            arguments
                options.species = string(missing)
                options.type = string(missing)
                options.CF0 = 20.
                options.l_factor = 3.8
                options.length = 20
                options.xs = string(missing)
                options.cfs = string(missing)
                options.rho = 1.000
                options.Ap = string(missing)
                options.bp = string(missing)
                options.Bu = string(missing)
                options.gain_const = string(missing)
                options.peak_magndb = string(missing)
                options.Bpeak = string(missing)
                options.phiaccum = string(missing)
                options.Nbeta = string(missing)
                options.Nf = string(missing)
                options.Qerb = string(missing)
                options.ERBbeta = string(missing)
                options.ERBf = string(missing)
                options.Qn = string(missing)
                options.Qn2 = string(missing)
                options.BWndBbeta = string(missing)
                options.BWndBf = string(missing)
                options.BWn2dBbeta = string(missing)
                options.BWn2dBf = string(missing)
                options.Sbeta = string(missing)
                options.Sf = string(missing)
                options.n = 10
                options.n2 = 3
                options.betas = string(missing)
                options.freqs = string(missing)
                options.aAp = string(missing)
                options.bAp = string(missing)
                options.aBu = string(missing)
                options.bBu = string(missing)
                options.pycochlea = string(missing)
            end
            if ~ismissing(options.pycochlea)
                pyCo = options.pycochlea;
            elseif ~ismissing(options.aAp)
                pyCo = py.Cochlea.Cochlea.five_param(type=options.type, ...
                    aAp=options.aAp, ...
                    bAp=options.bAp, ...
                    bp=options.bp, ...
                    aBu=options.aBu, ...
                    bBu=options.bBu, ...
                    gain_const=options.gain_const, ...
                    peak_magndb=options.peak_magndb, ...
                    CF0=options.CF0, ...
                    l_factor=options.l_factor, ...
                    length=options.length, ...
                    xs=options.xs, ...
                    cfs=options.cfs, ...
                    rho=options.rho, ...
                    betas=options.betas, ...
                    freqs=options.freqs);
            else
                pyCo = py.Cochlea.Cochlea(species=options.species, ...
                    type=options.type, ...
                    CF0=options.CF0, ...
                    l_factor=options.l_factor, ...
                    length=options.length, ...
                    xs=options.xs, ...
                    cfs=options.cfs, ...
                    rho=options.rho, ...
                    Ap=options.Ap, ...
                    bp=options.bp, ...
                    Bu=options.Bu, ...
                    gain_const=options.gain_const, ...
                    peak_magndb=options.peak_magndb, ...
                    Bpeak=options.Bpeak, ...
                    phiaccum=options.phiaccum, ...
                    Nbeta=options.Nbeta, ...
                    Nf=options.Nf, ...
                    Qerb=options.Qerb, ...
                    ERBbeta=options.ERBbeta, ...
                    ERBf=options.ERBf, ...
                    Qn=options.Qn, ...
                    Qn2=options.Qn2, ...
                    BWndBbeta=options.BWndBbeta, ...
                    BWndBf=options.BWndBf, ...
                    BWn2dBbeta=options.BWn2dBbeta, ...
                    BWn2dBf=options.BWn2dBf, ...
                    Sbeta=options.Sbeta, ...
                    Sf=options.Sf, ...
                    n=options.n, ...
                    n2=options.n2, ...
                    betas=options.betas, ...
                    freqs=options.freqs);
            end
            obj@FilterBank(pyfilterbank=pyCo)
            obj.PyCochlea = pyCo
        end
        function fil = filter_at_location(obj, x_coord, options)
            arguments
                obj
                x_coord
                options.gain_const = string(missing)
                options.peak_magndb = string(missing)
                options.type = 'P'
            end
            fil = Filter(pyfilter=obj.PyFilterBank.filter_at_location(x_coord, ...
                gain_const=options.gain_const, peak_magndb=options.peak_magndb, type=options.type));
        end
        function [betas, reals, imags, magns, phases] = plot_wavenumber(obj, options)
            arguments
                obj
                options.betas = string(missing)
                options.custom_title = 'Wavenumber (k)'
                options.show = true
                options.phase_in_rad = true
            end
            all_output = obj.PyFilterBank.plot_wavenumber(betas=options.betas, show=false, phase_in_rad=options.phase_in_rad);
            betas = double(all_output{1});
            reals = double(all_output{2});
            imags = double(all_output{3});
            magns = double(all_output{4});
            phases = double(all_output{5});
            
            if options.show
                plot(betas, reals)
                hold on
                plot(betas, imags)
                title(options.custom_title)
                xlabel('Normalized frequencies (Hz)')
                ylabel('Wavenumber (1/mm)')
                xline(obj.PyCochlea.bp_apex, ':k')
                yline(0, ':k')
            end
        end
        function [betas, reals, imags] = plot_impedance(obj, options)
            arguments
                obj
                options.betas = string(missing)
                options.custom_title = 'Normalized impedance (Z_{norm})'
                options.show = true
                options.phase_in_rad = true
            end
            all_output = obj.PyFilterBank.plot_wavenumber(betas=options.betas, show=false, phase_in_rad=options.phase_in_rad);
            betas = double(all_output{1});
            reals = double(all_output{2});
            imags = double(all_output{3});
            
            if options.show
                plot(betas, reals)
                hold on
                plot(betas, imags)
                title(options.custom_title)
                xlabel('Normalized frequencies (Hz)')
                ylabel('Normalized impedance (Ω/???)')
                xline(obj.PyCochlea.bp_apex, ':k')
                yline(0, ':k')
            end
        end
        function sigs = signal_response_heatmap(obj, signal, options)
            arguments
                obj
                signal
                options.len_xs = int32(20)
                options.custom_title = 'Signal Heatmap'
                options.show = true
            end
            % there must be a shorter way to do this
            data = obj.PyFilterBank.signal_response_heatmap(signal.PySignal, len_xs=int32(options.len_xs), show=false);
            signal.plot()
            int32(options.len_xs)
            data
            data = cellfun(@double, cell(data), UniformOutput=false);
            sigs = [];
            for row = data
                sigs = cat(1, sigs, cell2mat(row));
            end
            image(sigs);
        end
    end
end