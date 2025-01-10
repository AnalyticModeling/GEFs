classdef FilterBank
    properties
        PyFilterBank
        numFilters
        filters
    end
    methods
        function obj = FilterBank(options)
            arguments
                options.topology = string(missing)
                options.filters = string(missing)
                options.type = 'P'
                options.ir = string(missing)
                options.tf = string(missing)
                options.coeffs = string(missing)
                options.roots = string(missing)
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
                options.peak_freq = string(missing)
                options.pyfilterbank = string(missing)
            end
            if ismissing(options.pyfilterbank)
                pyFb = py.FilterBank.FilterBank(topology=options.topology, ...
                    filters=options.filters, ...
                    type=options.type, ...
                    ir=options.ir, ...
                    tf=options.tf, ...
                    coeffs=options.coeffs, ...
                    roots=options.roots, ...
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
                    peak_freq=options.peak_freq);
            else
                pyFb = options.pyfilterbank;
            end
            obj.PyFilterBank = pyFb;
            obj.numFilters = length(pyFb.filters);
            obj.filters = cellfun(@(x) Filter(pyfilter=x), cell(pyFb.filters));
        end
        function len = length(obj)
            len = obj.numFilters;
        end
        function fil = get_filter_from_uid(obj, uid)
            fil = obj.PyFilterBank.get_filter_from_uid(uid);
        end
        function sourceId = get_source_uid(obj, uid)
            sourceId = obj.PyFilterBank.get_source_uid(uid);
        end
        function childIds = get_uids_fed_into(obj, uid)
            childIds = obj.PyFilterBank.get_uids_fed_into(uid);
        end
        function add(obj, filter, opts)
            arguments
                obj
                filter
                opts.source = string(missing)
                opts.source_uid = -1
            end
            obj.PyFilterBank.add(filter, source=opts.source, source_uid=opts.source_uid)
        end
        function outsigs = process_signal(obj, signal, method)
            if nargin == 2
                method = 'default';
            end
            outsigs = OutputSignals(pyoutputsignals=obj.PyFilterBank.process_signal(signal.PySignal, method));
        end
        function allFilterData = bode_plot(obj, options)
            arguments
                obj
                options.freqs = string(missing)
                options.custom_title = 'Bode plot'
                options.show = true
            end
            allFilterData = obj.PyFilterBank.bode_plot(freqs=options.freqs, show=false);
            if options.show
                tiledlayout(2, 1, TileSpacing='tight')

                nexttile
                for i = 1:obj.numFilters
                    freqs = allFilterData{i}{1};
                    magns = allFilterData{i}{2};
                    uid = allFilterData{i}{4};
                    semilogx(freqs, magns)
                    hold on
                end
                % xlabel('Normalized frequency')
                ylabel('Magnitude (dB)')
                hold off

                nexttile
                for i = 1:obj.numFilters
                    freqs = allFilterData{i}{1};
                    phases = allFilterData{i}{3};
                    uid = allFilterData{i}{4};
                    semilogx(freqs, phases)
                    hold on
                end
                xlabel('Normalized frequency')
                ylabel('Phase (cycles)')
                hold off

                sgtitle(options.custom_title)
            end
        end
    end
end