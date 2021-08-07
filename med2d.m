function [y_final, f_final, kurtIter] = med2d(x,filterSize,termIter,termDelta,overlapMode,plotMode)
    %2D MINIMUM ENTROPY DECONVOLUTION WITH ADJUSTED CONVOLUTION FIX
    %  code by Geoff McDonald (glmcdona@gmail.com), February 2011
    %  Used in my MSc research at the University of Alberta, Advanced
    %  Control Systems Laboratory. This is a reference for my paper:
    %      McDonald, Geoff L., Qing Zhao, and Ming J. Zuo. "Maximum correlated Kurtosis deconvolution
    %      and application on gear tooth chip fault detection." Mechanical Systems and Signal Processing
    %      33 (2012): 237-255.
    %
    %  Updated in 2015 to include the convolution adjustment fix. Proposed
    %  in the second reference paper. This is important to fix MED from
    %  deconvolving the trivial solution at the convolution discontinuity.
    %  Using overlapMode of 'valid' uses the fix, while 'full' uses the
    %  original definition used by R.A. Wiggins (not recommended).
    %   
    %  Note to readers: This in general does not produce the optimal solution to
    %  the proposed deconvolution problem. You may want to consider the other implementation
    %  I've uploaded called 'omeda', which non-iteratively solves for the
    %  optimal solution to a similar problem based on the Minimum Entropy Deconvolution problem proposed
    %  by Carlos Cabrelli. Or if you are dealing with periodic impulses
    %  I recommend looking at 'momeda' which solves the optimal solution
    %  for deconvolving periodic impulses as a deconvolution goal.
    %
    %
    % med2d(x,filterSize,termIter,termDelta,overlapMode,plotMode)
    %
    % Algorithm Reference:
    %   Minimum Entropy Deconvolution:
    %    R.A. Wiggins, Minimum Entropy Deconvolution, Geoexploration, vol.
    %    16, Elsevier Scientific Publishing, Amsterdam, 1978. pp. 2135.
    %
    %   Convolution Adjustment Derivation:
    %    G.L. McDonald, Qing Zhao, Multipoint Optimal Minimum Entropy Deconvolution and Convolution
    %    Fix: Application to Vibration Fault Detection, unpublished
    %
    % Inputs:
    %    x: 
    %       Signal to perform Minimum Entropy Deconvolution on. If a single
    %       column/row of data is specified, a 1d filter is designed to
    %       minimize the entropy of the resulting signale. If a 2d data 
    %       matrix is specified, a single 1d filter will be designed to
    %       minimize the averaged entropy of each column of the filtered
    %       data.
    % 
    %    filterSize:
    %       This is the length of the finite inpulse filter filter to 
    %       design. Using a value of around 30 is appropriate depending on
    %       the data. Investigate the performance difference using
    %       different values.
    %
    %    termIter: (OPTIONAL)
    %       This is the termination number of iterations. If the 
    %       the number of iterations exceeds this number, the MED process
    %       will complete. Specify [] to use default value of 30.
    %
    %    termDelta: (OPTIONAL)
    %       This is the termination condition. If the change in kurtosis
    %       between iterations is below this threshold, the iterative
    %       process will terminate. Specify [] to use the default value
    %       of 0.01. You can specify a value of 0 to only terminate on
    %       the termIter condition, ie. execute an exact number of
    %       iterations.
    %
    %   overlapMode: (OPTIONAL)
    %       You should always use 'valid' for this parameter to include the
    %       convolution fix that corrects MED erroneously deconvolving 
    %       spurious impulses. See algorithm reference for Minimum Entropy
    %       Deconvolution Adjusted Convolution (MEDA) for details on the convolution
    %       adjustment. You can use 'full' if you want to reproduce the original
    %       MED method by R.A Wiggins, but it is not recommended for the above reason.
    % 
    %    plotMode:
    %       If this value is > 0, plots will be generated of the iterative
    %       performance and of the resulting signal.
    %
    % Outputs:
    %    y_final:
    %       The input signal(s) x, filtered by the resulting MED filter.
    %       This is obtained simply as: y_final = filter(f_final,1,x);
    %
    %    f_final:
    %       The final 1d MED filter in finite impulse response format.
    %
    %    kurtIter:
    %       Kurtosis according to MED iteration. kurtIter(end) is the
    %       final kurtosis, ie. the summed kurtosis of each y_final
    %       column of y_final. sum(kurtosis(each column of y_final))
    % 
    % Example:
    % % -------- 1d deconvolution example ------
    % n = 0:999;
    % x = [sin(n/30) + 0.2*(mod(n,21)==0)];
    % % Run for 100 iterations, 30 sample FIR filter
    % [y_final f_final kurt] = med2d(x',30,100,[],'valid',1);
    % 
    % % -------- 2d deconvolution example ------
    % % This will mostly extract the impulse-like
    % % disturbances caused by 0.2*(mod(n,21)==0)
    % % and plot the result.
    % n = 0:999;
    % x = [sin(n/30) + 0.2*(mod(n,21)==0);
    %   sin(n/13) + 0.2*(mod(n,21)==0)];
    % % Run for 100 iterations, 30 sample FIR filter
    % [y_final f_final kurt] = med2d(x',30,100,[],'valid',1);
    % 
    % 
    % Note:
    %   The solution is not the optimal solution to the
    %   entropy minimizataion problem, the solution is just a local
    %   minimum of the entropy and therefore a good pick.
    % Assign default values for inputs
    if( isempty(filterSize) )
        filterSize = 30;
    end
    if( isempty(termIter) )
        termIter = 100;
    end
    if( isempty(termDelta) )
        termDelta = -1000;
    end
    if( isempty(plotMode) )
        plotMode = 0;
    end
    
    % Validate the inputs
    if( sum( size(x) > 1 ) > 2 )
        error('MED:InvalidInput', 'Input signal x must be of either 2d or 1d.')
    %elseif( sum(size(termDelta) > 1) ~= 0 || termDelta < 0 )
    %    error('MED:InvalidInput', 'Input argument termDelta must be a positive scalar, or zero.')
    elseif( sum(size(termIter) > 1) ~= 0 || mod(termIter, 1) ~= 0 || termIter <= 0 )
        error('MED:InvalidInput', 'Input argument termIter must be a positive integer scalar.')
    elseif(  sum(size(plotMode) > 1) ~= 0 )
        error('MED:InvalidInput', 'Input argument plotMode must be a scalar.')
    elseif( sum(size(filterSize) > 1) ~= 0 || filterSize <= 0 || mod(filterSize, 1) ~= 0 )
        error('MED:InvalidInput', 'Input argument filterSize must be a positive integer scalar.')
    end
    
    % Validate the inputs
    if( strcmp(overlapMode,'full') == 1 )
        overlap_full = 1; % not recommended
        warning('MED WARNING: Using the full option often erroneously deconvolves an impulse near the start of the signal. It is recommended to use the valid option instead')
    elseif( strcmp(overlapMode,'valid') == 1 )
        overlap_full = 0;
    else
        error('MEDD:NoOverlapMode', 'overlapMode argument must be "valid" or "full".')
    end
    % If the data is 1d, lets make it a column vector
    if( sum(size(x)>1) == 1 )
        x = x(:); % A column vector
    end
    L = filterSize;
    
    %%% First we need to calculate the autocorrelation matrix R
    N = length(x);
    Xm0 = zeros(L,N+L-1); % y = f*x where x is padded
    
    for( l =1:L )
        if( l == 1 )
            Xm0(l,1:N) = x(1:N);
        else
            Xm0(l,2:end) = Xm0(l-1, 1:end-1);
        end
    end
    
    if( ~overlap_full )     % "valid" region only
        Xm0 = Xm0(:,L:N-1); % y = f*x where only valid x is used
                            % y = Xm0'x to get valid output signal
    end
    
    average_autocorr = Xm0*Xm0';
    
    % Generate the toeplitz matrix R and it's inverse
    R_inv = pinv(average_autocorr);
    
    
    % Initialize matrix sizes
    f = zeros(L,1);
    y = zeros(size(x,1),size(x,2));
    b = zeros(L,1);
    kurtIter = [];
    
    % Assume initial filter as a difference filter
    f(2) = 1;
    
    % Iteratively adjust the filter to minimize entropy
    n = 1;
    while n == 1 || ( n < termIter && ( (kurt(filter(f,1,x)) - kurtIter(n-1)) > termDelta ) )
        % Compute output signal
        y = Xm0'*f;
        
        % Calculate the kurtosis
        kurtIter(n) = kurt(y); %#ok<AGROW>
        
        yc = y.^3;
        weightedCrossCorr = Xm0*yc;
        
        % Now we have new filter coefficients calculted as:
        % f = A^-1 * g
        f = R_inv*weightedCrossCorr;
        f = f/sqrt(sum(f.^2)); % Normalize the filter result
        % Next iteration
        n = n + 1;
    end
    
    % Update the final result
    f_final = f;
    y_final = Xm0'*f_final;
    kurtIter(n) = kurt(y_final);
    
    % Plot the results
    if( plotMode > 0 )
        
        figure;
        subplot(2,1,1)
        plot(x)
        title('Input Signal(s)')
        ylabel('Value')
        xlabel('Sample Number')
        axis tight
        
        subplot(2,1,2)
        plot(y_final)
        title('Signal(s) Filtered by MED')
        ylabel('Value')
        xlabel('Sample Number')
        axis tight
        
        figure;
        stem(f_final)
        xlabel('Sample Number')
        ylabel('Value')
        title('Final Filter, Finite Impulse Response')
        figure;
        plot(kurtIter);
        xlabel('MED Algorithm Iteration')
        ylabel('Sum of Kurtosis for Filtered Signal(s)')
        
        if( n == termIter )
            display('Terminated for iteration condition.')
        else
            display('Terminated for minimum change in kurtosis condition.')
        end
    end
end
function [result] = kurt(x)
    % This function simply calculates the summed kurtosis of the input
    % signal, x.
    result = mean( (sum((x-ones(size(x,1),1)*mean(x)).^4)/(size(x,1)-2))./(std(x).^4) );
    %result = sum(x.^4)/(sum(x.^2)^2);
end