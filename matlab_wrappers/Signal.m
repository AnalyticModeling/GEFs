classdef Signal
    properties
        PySignal
        length
    end
    methods
        function obj = Signal(options)
            arguments
                options.mode = 't'
                options.data = ones(1, 9)
                options.fs = 1
                options.evenlen = true
                options.func = string(missing)
                options.num_samples = 9
                options.f_init = string(missing)
                options.w_init = string(missing)
                options.f_final = string(missing)
                options.w_final = string(missing)
                options.freq_func = string(missing)
                options.freqs = string(missing)
                options.init_phase = string(missing)
                options.pysignal = string(missing)
            end
            if ismissing(options.pysignal)
                if ~ismissing(options.func)
                    pySig = py.Signal.Signal.from_function(mode=options.mode, func=options.func, fs=options.fs, evenlen=options.evenlen);
                elseif ~ismissing(options.f_init) || ~ismissing(options.w_init)
                    pySig = py.Signal.Signal.linear_chirp(f_init=options.f_init, w_init=options.w_init, f_final=options.f_final, w_final=options.w_final, fs=options.fs, num_samples=options.num_samples);
                elseif ~ismissing(options.freq_func)
                    pySig = py.Signal.Signal.from_instantaneous_frequency(freq_func=options.freq_func, freqs=options.freqs, init_phase=options.init_phase, fs=options.fs, num_samples=options.num_samples);
                else
                    pySig = py.Signal.Signal(mode=options.mode, data=options.data, fs=options.fs);
                end
            else
                pySig = options.pysignal;
            end
            obj.PySignal = pySig;
            obj.length = pySig.length;
        end
        function r = plus(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1+s2', 'res', s1=s1, s2=s2));
        end
        function r = minus(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1-s2', 'res', s1=s1, s2=s2));
        end
        function r = uminus(a)
            r = Signal(pysignal=pyrun('res=-s1', 'res', s1=a.PySignal));
        end
        function r = uplus(a)
            r = Signal(pysignal=pyrun('res=+s1', 'res', s1=a.PySignal));
        end
        function r = times(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1*s2', 'res', s1=s1, s2=s2));
        end
        function r = rdivide(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1/s2', 'res', s1=s1, s2=s2));
        end
        function r = ldivide(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s2/s1', 'res', s1=s1, s2=s2));
        end
        function r = power(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1**s2', 'res', s1=s1, s2=s2));
        end
        function r = floorDiv(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1//s2', 'res', s1=s1, s2=s2));
        end
        function r = mod(a, b)
            if isa(a, 'Signal'); s1 = a.PySignal; else; s1 = a; end
            if isa(b, 'Signal'); s2 = b.PySignal; else; s2 = b; end
            r = Signal(pysignal=pyrun('res=s1%s2', 'res', s1=s1, s2=s2));
        end
        function r = abs(a)
            r = Signal(pysignal=pyrun('res=abs(s1)', 'res', s1=a.PySignal));
        end
        function val = at_time(obj, t, tolerance)
            if nargin == 2
                tolerance = 1e-10;
            end
            val = obj.PySignal.at_time(t, tolerance);
        end
        function dataseries = get_data(obj, mode)
            if nargin == 1
                mode = 't';
            end
            dataseries = double(obj.PySignal.get_data(mode));
        end
        function newsig = resample(obj, new_sample_rate, end_time)
            if nargin == 2
                newsig = Signal(pysignal=obj.PySignal.resample(new_sample_rate));
            else
                newsig = Signal(pysignal=obj.PySignal.resample(new_sample_rate, end_time));
            end
        end
        function [upper, lower] = envelope_analytic(obj)
            all_output = obj.PySignal.envelope_analytic();
            upper = double(all_output{1});
            lower = double(all_output{2});
        end
        function instphi = instantaneous_phase(obj)
            instphi = double(obj.PySignal.instantaneous_phase());
        end
        function instf = instantaneous_freq(obj)
            instf = double(obj.PySignal.instantaneous_freq());
        end
        function specH = spectral_entropy(obj)
            specH = double(obj.PySignal.spectral_entropy());
        end
        function moveH = moving_spectral_entropy(obj, window_len)
            if nargin == 1
                window_len = 9;
            end
            moveH = double(obj.PySignal.moving_spectral_entropy(window_len=window_len));
        end
        function [sfft, bounds] = spectrogram(obj, options)
            arguments
                obj
                options.win = string(missing)
                options.hop = int32(1)
                options.mfft = int32(200)
                options.custom_title = 'Spectrogram'
                options.show = true
            end
            output = obj.PySignal.spectrogram(win=options.win, hop=int32(options.hop), mfft=int32(options.mfft), ...
                custom_title=options.custom_title, show=options.show);
            output
            sfft = 0;
            bounds = 0;
        end
        function r = crosscorrelate(a, b)
            r = double(pyrun('res=s1.crosscorrelate(s2)', 'res', s1=a.PySignal, s2=b.PySignal));
        end
        function corr = autocorrelate(obj, options)
            arguments
                obj
                options.custom_title = 'Autocorrelation'
                options.show = true
            end
            corr = obj.PySignal.autocorrelate(show=false);
            % incomplete
        end
        function plot(obj)
            plot(obj.PySignal.timestamps, obj.PySignal.mode_t)
        end
    end
end