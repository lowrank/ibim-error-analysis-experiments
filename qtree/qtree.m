
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QTREE
%
% Generate kd-tree structure for implicit boundary integral method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef qtree < handle
    %QTREE kd-tree

    properties (Access =public)

        dim

        dict
        maxId
        root
        nLevelMax
        h

        shift
        rotation

        distFunc
        EPS

        validList
    end

    properties (Access = private)


        center_diffs

    end

    methods
        function obj = qtree(d_, h_, EPS_,randShift, randRotation )
            %QTREE Construct kd-tree
            if nargin < 5
                randRotation = false;
            end

            if nargin < 4
                randShift = false;
            end

            obj.dim    = d_;
            obj.h      = h_;
            obj.EPS    = EPS_;
            obj.root   = 1;
            obj.nLevelMax = 1 + ceil(-log(h_) / log(2)); % height of tree
            
            
            obj.maxId = obj.root;
            obj.validList = [];

            if randShift
                obj.shift = ( rand(obj.dim, 1) - 0.5 ) * h_;
            else
                obj.shift =  zeros(obj.dim, 1) ;
            end

            if randRotation
                obj.rotation = RandUGroup(obj.dim, 'GType', 'SO');
            else 
                obj.rotation = eye(obj.dim);
            end

            % prepare private properties
            obj.center_diffs = zeros(obj.dim, 2^obj.dim);
            for i = 1: 2^obj.dim
                for dimIndex = 1:obj.dim
                     obj.center_diffs(dimIndex, i) = (  bitand( bitshift(i-1,  -dimIndex+1), 1) - 0.5  );
                end
            end
           
        end

        function ret = check(obj, n)
            % The distance function is necessary.

            d_ =  obj.distFunc(n.center);
            if sqrt(obj.dim) * n.radius > abs(d_)
                ret = true;
            elseif abs(d_) - sqrt(obj.dim) * n.radius < obj.EPS
                ret = true;
            else
                ret = false;
            end
        end

        function obj = populate(obj)
            %POPULATE Generate kd-tree
            obj.root = 1;    % first node is root.
            obj.dict{obj.root} = qnode(obj.dim, 1, 1);   % level = 1, index = 1
            obj.dict{obj.root}.center = obj.shift;
            obj.dict{obj.root}.radius = 2^(obj.nLevelMax - 1) * obj.h;
            obj.assignNodes(obj.root);
            
        end

        function obj = assignNodes(obj, id_)
            if obj.dict{id_}.nLevel == obj.nLevelMax+1
                obj.dict{id_}.isLeaf = true;
                dist_ =  obj.distFunc(obj.dict{id_}.center);

                if abs (dist_) < obj.EPS  
                    obj.dict{id_}.isValid = true;
                    obj.dict{id_}.dist = dist_; 
                    obj.validList(end + 1) = id_;
                end
            else % decide if needs division
                
                if obj.check(obj.dict{id_}) 

                    % division
                    for i = 1:(2^obj.dim)
                        obj.maxId = obj.maxId + 1;
                        obj.dict{id_}.child(i) = obj.maxId;

                        obj.dict{ obj.maxId } =  qnode(obj.dim, i, obj.dict{id_}.nLevel + 1);
                        obj.dict{ obj.maxId }.parent = id_;

                        % add properties: center, radius

                        obj.dict{ obj.maxId }.center =...
                                obj.dict{ id_ }.center + ...
                                obj.rotation * obj.dict{ id_ }.radius *...
                                obj.center_diffs(:, i);

                        obj.dict{ obj.maxId }.radius = obj.dict{ id_ }.radius / 2;

                    end
                    % recursive calls
                    for i = 1:(2^obj.dim)
                        obj = assignNodes(obj, obj.dict{id_}.child(i));
                    end                    
                end
            end
        end
    end


end