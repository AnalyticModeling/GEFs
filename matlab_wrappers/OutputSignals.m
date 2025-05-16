classdef OutputSignals
    properties
        PyOutputSignals
        numSignals
        inSignal
        outSignals
    end
    methods
        function obj = OutputSignals(options)
            arguments
                options.all_signals = string(missing)
                options.graph = string(missing)
                options.pyoutputsignals = string(missing)
            end
            if ismissing(options.pyoutputsignals)
                pyOs = py.OutputSignals.OutputSignals(options.all_signals, options.graph);
            else
                pyOs = options.pyoutputsignals;
            end
            obj.PyOutputSignals = pyOs;
            obj.numSignals = length(pyOs.outsignals);
            obj.inSignal = Signal(pysignal=pyOs.insignal);
            obj.outSignals = cellfun(@(x) Signal(pysignal=x), cell(pyOs.outsignals));
        end
        function len = length(obj)
            len = obj.numSignals;
        end
        function sig = get_signal_from_uid(obj, uid)
            sig = obj.PyOutputSignals.get_signal_from_uid(uid);
        end
        function sourceId = get_source_uid(obj, uid)
            sourceId = obj.PyOutputSignals.get_source_uid(uid);
        end
        function correlogram(obj, custom_title)
            if nargin == 1
                custom_title = 'Correlogram';
            end
            n = obj.numSignals-1;
            sigLen = int32(obj.PyOutputSignals.signal_length);
            obj.PyOutputSignals
            fullxaxis = double(-sigLen+1:sigLen-1) ...
                ./ double(obj.PyOutputSignals.signal_fs);
            tiledlayout(n, n, TileSpacing='tight')
            for y = 2:n+1
                for x = 1:n
                    if y >= x+1
                        nexttile(n*(y-2)+x)
                        corr = crosscorrelate(obj.outSignals(y), obj.outSignals(x));
                        plot(fullxaxis, corr)
                        xticklabels([])
                    end
                end
            end
            for i = 1:n
                nexttile(1+(i-1)*n)
                ylabel(sprintf('Signal %u', obj.outSignals(i+1).PySignal.uid))
                nexttile((n-1)*n + i)
                xlabel(sprintf('Signal %u', obj.outSignals(i).PySignal.uid))
            end
            sgtitle(custom_title)
        end
        function autocorrelates(obj, custom_title)
            if nargin == 1
                custom_title = 'Autocorrelates';
            end
            tiledlayout(obj.numSignals, 1, TileSpacing='tight')
            xlim = int32(obj.PyOutputSignals.signal_length);
            for i = 1:obj.numSignals
                nexttile
                plot(0:xlim-1, autocorrelate(obj.outSignals(i)))
                ylabel(sprintf('Signal %u', obj.outSignals(i).PySignal.uid))
            end
            sgtitle(custom_title)
        end
        function correlate_with(obj, signal, custom_title)
            if nargin == 2
                custom_title = 'Correlations';
            end
            tiledlayout(obj.numSignals, 1, TileSpacing='tight')
            xlim = int32(obj.PyOutputSignals.signal_length);
            for i = 1:obj.numSignals
                nexttile
                size(crosscorrelate(signal, obj.outSignals(i)))
                plot(-xlim+1:xlim-1, crosscorrelate(signal, obj.outSignals(i)))
                ylabel(sprintf('Signal %u', obj.outSignals(i).PySignal.uid))
            end
            sgtitle(custom_title)
        end
    end
end