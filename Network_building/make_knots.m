function [ cntrl_pts ] = make_knots( nlags, number_knots )
% MAKE_KNOTS creates a vector of nonuniform numbers that can be used to
%            work on knot placement, 'number_knots' includes last at 'nlags'
%            and excludes first at '0'

% cntrl_pts = [0 nlags];
% indx = nlags;
% 
%     for k = 1:number_knots-2
%         indx = round(.5*(cntrl_pts(1)+indx));
%         cntrl_pts = [cntrl_pts indx];
%     end
% 
% cntrl_pts = sort(unique(cntrl_pts));

closeness = 2;
cntrl_pts = [0 1:closeness-1 nlags];

remaining_number_knots = number_knots - length(cntrl_pts);
distance = nlags - closeness - 2;
linear_spacing = floor(distance/remaining_number_knots);

cntrl_pts = [cntrl_pts (closeness:linear_spacing:nlags)];

cntrl_pts = sort(unique(cntrl_pts),'ascend');

if cntrl_pts(end)-cntrl_pts(end-1) < 10
    cntrl_pts(end-1) = [];
end

end

