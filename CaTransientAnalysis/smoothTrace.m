function [traceSmoothed, traceSmoothedDetrended] = smoothTrace(tracePixel, filterWidth, varargin)
            % This function smooths (~denoises) a trace using Savitzky-Golay filter of
            % a given width. It may also remove drift/trend in data using
            % polynomial fitting.
            %
            % IN:
            % tracePixel - the trace to be smoothed.
            %
            % filterWidth - the width of the Savitzky-Golay filter (should
            % be an odd number - if an even one is provided, one is added
            % to it).
            %
            % varargin{1} - the optional parameter may contain the order of the polynomial
            % used to detrend the signal (useful for removing drift etc.).
            % Unlike standard methods which subtract the polynomial making
            % the signal approximately zero-centered on the y-axis, here we
            % re-add the average of the original trace to the zero-centered
            % trace, so that the information on signal baseline is not
            % lost.
            %
            % OUT:
            % traceSmoothed - tracePixel after smoothing.
            %
            % traceSmoothedDetrended - tracePixel after smoothing and
            % baseline subtraction.
            
            tracePixel = tracePixel(:); % making sure this is a column vector
            %% detrending
            toDetrend = ~isempty(varargin);
            
            if toDetrend
                detrendOrder = varargin{1};
                assert(isnumeric(detrendOrder), 'Order of baseline subtraction must be numeric');
                assert(detrendOrder>0, 'Order of baseline subtraction must be larger than 0');
                
                % below, we fit a polynomial to the signal, which acts as a
                % baseline. Mu shenanigans serve to improve the numerical
                % stability (see documentation of polyfit for mu).
                [p,s,mu] = polyfit([1:length(tracePixel)]', tracePixel, detrendOrder); % mean subtraction is just so that the problem isn't ill-conditioned
                polyBaseline = polyval(p, [1:length(tracePixel)]', [], mu);
                traceDeTrend = tracePixel - polyBaseline + mean(tracePixel); %we re-add mean of tracePixel so that we don't lose information on total signal baseline
            end
            
            if (mod(filterWidth, 2) == 0)
                filterWidth = filterWidth + 1;
            end
            
            if (filterWidth > 4) % if filter width is <= 4, it's not possible to carry out filtering because the filter width must be more than filter order. In this way, we can also forbid filtering (by setting smoothingParameter to 1)
                traceSmoothed = sgolayfilt(tracePixel,4,filterWidth);
                traceSmoothedDetrended = traceSmoothed; % just so that this is defined as an output parameter if no detrending is used.
                if toDetrend
                    traceSmoothedDetrended = sgolayfilt(traceDeTrend,4,filterWidth);
                end
            else
                traceSmoothed = tracePixel;
                traceSmoothedDetrended = traceSmoothed;
                if toDetrend
                    traceSmoothedDetrended = traceDeTrend;
                end
            end
            
            
        end
