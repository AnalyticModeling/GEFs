classdef Filter
    properties
        PyFilter
    end
    methods
        function obj = Filter(options)
            arguments
                options.ir = string(missing)
                options.tf = string(missing)
                options.coeffs = string(missing)
                options.roots = string(missing)
                options.type = 'P'
                options.Ap = string(missing)
                options.bp = string(missing)
                options.Bu = string(missing)
                options.gain_const = string(missing)
                options.peak_magndb = 0
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
                options.cf = string(missing)
                options.pyfilter = string(missing)
            end
            if ismissing(options.pyfilter)
                pyFil = py.Filter.Filter(ir=options.ir, ...
                    tf=options.tf, ...
                    coeffs=options.coeffs, ...
                    roots=options.roots, ...
                    type=options.type, ...
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
                    freqs=options.freqs, ...
                    cf=options.cf);
            % multiband?
            else
                pyFil = options.pyfilter;
            end
            obj.PyFilter = pyFil;
        end
        function cs = get_computed_chars(obj)
            cs = struct(obj.PyFilter.get_computed_chars());
        end
        function cs = get_computed_unnormalized_chars(obj)
            cs = struct(obj.PyFilter.get_computed_unnormalized_chars());
        end
        function cs = get_orig_chars(obj)
            cs = struct(obj.PyFilter.get_orig_chars());
        end
        function cs = get_params(obj)
            cs = struct(obj.PyFilter.get_params());
        end
        function cs = solve(obj, sig, options)
            arguments
                obj
                sig
                options.method = string(missing)
                options.fs = 100
            end
            cs = Signal(pysignal=obj.PyFilter.solve(sig.PySignal, method=options.method, fs=options.fs));
        end
        function [freqs, magns, phases] = bode_plot(obj, options)
            arguments
                obj
                options.freqs = string(missing)
                options.custom_title = 'Bode plot'
                options.show = true
            end
            all_output = obj.PyFilter.bode_plot(freqs=options.freqs, show=false);
            freqs = double(all_output{1});
            magns = double(all_output{2});
            phases = double(all_output{3});
            if options.show
                tiledlayout(2, 1, TileSpacing='tight')

                nexttile
                semilogx(freqs, magns)
                xlabel('Normalized frequency')
                ylabel('Magnitude (dB)')

                nexttile
                semilogx(freqs, phases)
                xlabel('Normalized frequency')
                ylabel('Phase (cycles)')

                sgtitle(options.custom_title)
            end
        end
        function [freqs, reals, imags] = frequency_real_imag_plot(obj, options)
            arguments
                obj
                options.num_samples = int32(10000)
                options.freqs = string(missing)
                options.peak_magn = 1
                options.custom_title = 'Frequency response plot'
                options.show = true
            end
            all_output = obj.PyFilter.frequency_real_imag_plot(num_samples=int32(options.num_samples), ...
                    freqs=options.freqs, ...
                    peak_magn=options.peak_magn, ...
                    show=false);
            freqs = double(all_output{1});
            reals = double(all_output{2});
            imags = double(all_output{3});
            if options.show
                tiledlayout(2, 1, TileSpacing='tight')

                nexttile
                semilogx(freqs, reals)
                xlabel('Normalized frequency')
                ylabel('Re(z)')

                nexttile
                semilogx(freqs, imags)
                xlabel('Normalized frequency')
                ylabel('Im(z)')

                sgtitle(options.custom_title)
            end
        end
        function [freqs, magns, phases] = nichols_plot(obj, options)
            arguments
                obj
                options.num_samples = int32(10000)
                options.freqs = string(missing)
                options.custom_title = 'Nichols plot'
                options.show = true
            end
            all_output = obj.PyFilter.bode_plot(num_samples=int32(options.num_samples), ...
                    freqs=options.freqs, ...
                    show=false);
            freqs = double(all_output{1});
            magns = double(all_output{2});
            phases = double(all_output{3});
            if options.show
                plot(phases, magns)
                xlabel('Phase (rad)')
                ylabel('Magnitude (dB)')
                title(options.custom_title)
            end
        end
        function [freqs, reals, imags] = nyquist_plot(obj, options)
            arguments
                obj
                options.num_samples = int32(10000)
                options.freqs = string(missing)
                options.custom_title = 'Nyquist plot'
                options.show = true
            end
            all_output = obj.PyFilter.frequency_real_imag_plot(num_samples=int32(options.num_samples), ...
                    freqs=options.freqs, ...
                    show=false);
            freqs = double(all_output{1});
            reals = double(all_output{2});
            imags = double(all_output{3});
            if options.show
                plot(reals, imags)
                xlabel('Re(z)')
                ylabel('Im(z)')
                title(options.custom_title)
            end
        end
        function [timestamps, response] = impulse_response_plot(obj, options)
            arguments
                obj
                options.num_samples = int32(10000)
                options.sample_freq = 500
                options.custom_title = 'Impulse response'
                options.show = true
            end
            all_output = obj.PyFilter.impulse_response_plot(num_samples=int32(options.num_samples), ...
                    sample_freq=options.sample_freq, ...
                    show=false);
            timestamps = double(all_output{1});
            response = double(all_output{2});
            if options.show
                plot(timestamps, response)
                xlabel('Time (s)')
                title(options.custom_title)
            end
        end
        function [zeros, poles] = pzplot(obj, options)
            arguments
                obj
                options.custom_title = string(missing)
                options.show = true
            end
            all_output = obj.PyFilter.pole_zero_plot(show=false);
            zeros = double(all_output{1});
            poles = double(all_output{2});
            if options.show
                if ismissing(options.custom_title)
                    try
                        Bu = obj.get_params().Bu;
                        if abs(Bu-round(Bu))>1e-10
                            custom_title = 'Pole-zero plot (of base filter)';
                        else
                            custom_title = 'Pole-zero plot';
                        end
                    catch
                        custom_title = 'Pole-zero plot';
                    end
                else
                    custom_title = options.custom_title;
                end
                zeroreal = real(zeros);
                zeroimag = imag(zeros);
                scatter(zeroreal, zeroimag, [], 'blue', 'x', LineWidth=2)
                hold on
                polereal = real(poles);
                poleimag = imag(poles);
                scatter(polereal, poleimag, [], [1 0.5 0], 'o', LineWidth=2)
                xline(0)
                yline(0)
                xlabel('Re(z)')
                ylabel('Im(z)')
                axis equal
                title(custom_title)
            end
        end
        function [xaxis, Qns] = Qn_plot(obj, options)
            arguments
                obj
                options.max_n = int32(20)
                options.custom_title = 'Qn plot'
                options.show = true
            end
            all_output = obj.PyFilter.Qn_plot(max_n=int32(options.max_n), show=false);
            xaxis = double(all_output{1});
            Qns = double(all_output{2});
            if options.show
                plot(xaxis, Qns)
                xlabel('n (dB)')
                ylabel('Q_n')
                title(options.custom_title)
            end
        end
        function [errorableChars, errors] = characteristic_error(obj, options) % is xaxis a good name
            arguments
                obj
                options.custom_title = 'Estimated vs Desired Bar Chart'
                options.show = true
            end
            all_output = obj.PyFilter.characteristic_error(show=false);
            errorableChars = string(all_output{1});
            errors = double(all_output{2});
            if options.show
                bar(errorableChars, errors)
                title(options.custom_title)
            end
        end
    end
end