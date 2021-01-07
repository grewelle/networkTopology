function graph_mwb(loc, H, mwb,n, A, mobility)
    colormap = jet(); %256 by 3 array
    if max(mwb) >30
        max_mwb = max(mwb);
    else
        max_mwb = 30;
    end
    indices = ceil(mwb.*((64-1)./max_mwb) + 1)'; %scale colors by mwb
    colormap(indices,:);
    
    figure; 
    
    %plot lines between connected nodes
    for i = 1:n
        for j = 1:n
            if A(i,j) > 0
                line([loc(i,1) loc(j,1)],[loc(i,2) loc(j,2)],'color', 'k');
                %hold on;

            end
        end
    end
    maptitle = sprintf('Map of Mean Worm Burden at m = %d', mobility);
    title(maptitle);

    %plot circles representing each node. The radius is proportional to the
    %population size
    for node = 1:n
       % hold on;
        circles(loc(node,1),loc(node,2),log10(H(node))/2+10, 'Color', colormap(indices(node),:));
    end
  
    labels = cellstr( num2str([1:size(loc,1)]') );  %' # labels correspond to their order
    text(loc(:,1)+5, loc(:,2)+5, labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right')
    axis equal; axis off;
    %hold off;
    %gplot(A,loc)
end