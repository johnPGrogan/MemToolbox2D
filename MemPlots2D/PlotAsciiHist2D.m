% PlotAsciiHist2D makes a cute little ascii histogram of error data
%
%   PlotAsciiHist2D(data,n);
%
% n says how many bins to use. So PlotAsciiHist2D(data, 12) prints something
% like this:
%
%     X:___.-'-.___  Y:___.-'-.___
% 
% This 2D version is adapted from code in the MemToolbox package(Suchow, J. W.,
% Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). Modeling visual working 
% memory with the MemToolbox. Journal of Vision, 13(10):9, 1–8. 
% doi:10.1167/13.10.9. MemToolbox.org). 
% It was adapted by John P. Grogan, 2019.
%

function PlotAsciiHist2D(data,n)
    % Check input
    if nargin < 2
        n = 21;
    end

    if size(data,1) > 1 %if 2D data
        % Bin the data
        m = hist3(data', [n,n]);

        %take average over X and Y dimensions
        mXY(1,:) = max(m,[],2);
        mXY(2,:) = max(m,[],1);

        % Figure out highest bin
        maxBin = [max(mXY(1,:)),max(mXY(2,:))];
    
        % Build the histogram
        symbols = {'_', '.', '-', ''''};
        h = {'X:','Y:'};
        for j = 1:2
            for i = 1:n
                symbolToAdd = symbols{1+floor((mXY(j,i)/maxBin(j))*(length(symbols)-1))};
                h{j} = [h{j},symbolToAdd];
            end
            h{j} = [h{j}, ' '];
        end
        
    else
        error('Your data is not 2D')
    end

    % Display the histogram
    disp(cat(2,h{:}));
end

