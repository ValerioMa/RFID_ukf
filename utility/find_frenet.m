function [P_frenet] = find_frenet(path,robot)
%FIND_FRENET Summary of this function goes here
%   Detailed explanation goes here

deltas = path(1:2,:)-robot(1:2);
[~,idx]=min(deltas(1,:).*deltas(1,:)+deltas(2,:).*deltas(2,:));
P_frenet = path(:,idx);
end

