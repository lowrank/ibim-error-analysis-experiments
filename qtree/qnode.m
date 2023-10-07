
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QNODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef qnode < handle
    % QNODE Defines node class for kd-tree structure in 2d/3d.

    properties(Access = public)
        dim
        center
        radius

        parent
        child

        nLevel
        nIndex

        isLeaf
        isValid

        dist
    end

    methods
        function obj = qnode(d_, index_, level_)
            %QNODE Construct a node (initialized)
            obj.dim     = d_;
            obj.nLevel  = level_;
            obj.nIndex  = index_;
            obj.parent  = -1;                         % no parent initially
            obj.isValid = false;
            obj.isLeaf  = false;
            obj.child   = -1 * ones(2^obj.dim, 1);
            obj.center  = -inf * ones(obj.dim, 1);
            obj.radius  = inf;
            obj.dist    = -inf;              % dist from center to boundary
        end
    end
end