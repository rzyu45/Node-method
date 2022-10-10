classdef NMpipe
    % Pipes are actually queues in node method

    properties
        water_mass %queues of water mass
        T%queues of temperature
        m_buffer
        T_buffer
        parameter
        index
        num
        opt
        v % total mass flow in the pipe. unit: kg
        T_updated % bool if the inlet temperature has been updated
        direction % mass flow direction
    end

    methods
        function obj = NMpipe(parameter, index, num, opt, m0, T0)
            % Constructor function
            
            obj.parameter=parameter;
            obj.opt=opt;
            obj.num=num.pipe;
            obj.index=index;
            obj.water_mass=cell(obj.num,1);
            obj.T=cell(obj.num,1);
            obj.v=parameter.density.*parameter.area.*parameter.L;
            obj.T_updated=true;
            obj.direction=ones(obj.num,1);
            obj.m_buffer=cell(obj.num,1);
            obj.T_buffer=cell(obj.num,1);

            for pipe=obj.index.pipe.'
                obj.water_mass{pipe}=m0(pipe)*ones(1,ceil(obj.v(pipe)/m0(pipe)/opt.t));
                obj.T{pipe}=T0(pipe)*ones(size(obj.water_mass{pipe}));
                obj.m_buffer{pipe}=m0(pipe);
                obj.T_buffer{pipe}=T0(pipe);
            end

        end

        function obj=push(obj,m)
            if obj.T_updated==true
                for pipe=obj.index.pipe.'
                    if obj.direction(pipe)==1
                        obj.water_mass{pipe}=[m(pipe) obj.water_mass{pipe}];
                        % push a fake T because at this time the true T has not been calculated
                        obj.T{pipe}=[obj.T{pipe}(1) obj.T{pipe}];
                    else
                        obj.water_mass{pipe}=[obj.water_mass{pipe} m(pipe)];
                        % push a fake T because at this time the true T has not been calculated
                        obj.T{pipe}=[obj.T{pipe} obj.T{pipe}(end)];
                    end
                end
                obj.T_updated=false;
                obj=obj.get_buffer();
            else
                error("Inlet temperature not updated! Push new water mass prohibited!")
            end
        end

        function [obj, m , T]=pop(obj, pipe)
                if obj.direction(pipe)==1
                    m=obj.water_mass{pipe}(end);
                    obj.water_mass{pipe}(end)=[];
                    T=obj.T{pipe}(end);
                    obj.T{pipe}(end)=[];
                else
                    m=obj.water_mass{pipe}(1);
                    obj.water_mass{pipe}(1)=[];
                    T=obj.T{pipe}(1);
                    obj.T{pipe}(1)=[];
                end
        end
        
        function [m,T]=get_outlet(obj,pipe)
            if obj.direction(pipe)==1
                m=obj.water_mass{pipe}(end);
                T=obj.T{pipe}(end);
            else
                m=obj.water_mass{pipe}(1);
                T=obj.T{pipe}(1);
            end
        end

        function [m,T]=get_inlet(obj,pipe)
            if obj.direction(pipe)==1
                m=obj.water_mass{pipe}(1);
                T=obj.T{pipe}(1);
            else
                m=obj.water_mass{pipe}(end);
                T=obj.T{pipe}(end);
            end
        end

        function obj=get_buffer(obj)
            obj=obj.clean_buffer();
            for pipe=obj.index.pipe.'
                [m_outlet,T_outlet]=get_outlet(obj,pipe);
                [m_inlet,~]=get_inlet(obj,pipe);
                while sum(obj.water_mass{pipe}*obj.opt.t)-m_outlet*obj.opt.t>obj.v(pipe)+m_inlet*obj.opt.t
                    [obj,~,~]=obj.pop(pipe);
                    [m_outlet,T_outlet]=get_outlet(obj,pipe);
                end
                if sum(obj.water_mass{pipe}*obj.opt.t)-m_outlet*obj.opt.t<=obj.v(pipe)
                    obj=add_T_buffer(obj,T_outlet,pipe);
                    obj=add_m_buffer(obj,m_inlet,pipe);
                else
                    obj=add_T_buffer(obj,T_outlet,pipe);
                    obj=add_m_buffer(obj,(m_outlet*obj.opt.t-(sum(obj.water_mass{pipe}*obj.opt.t)-obj.v(pipe)-m_inlet*obj.opt.t))/obj.opt.t,pipe);
                    [obj,~,~]=obj.pop(pipe);
                    [m_outlet,T_outlet]=get_outlet(obj,pipe);
                    while sum(obj.water_mass{pipe}*obj.opt.t)-m_outlet*obj.opt.t>obj.v(pipe)
                        obj=add_T_buffer(obj,T_outlet,pipe);
                        obj=add_m_buffer(obj,m_outlet,pipe);
                        [obj,~,~]=obj.pop(pipe);
                        [m_outlet,T_outlet]=get_outlet(obj,pipe);
                    end
                    obj=add_T_buffer(obj,T_outlet,pipe);
                    obj=add_m_buffer(obj,(sum(obj.water_mass{pipe}*obj.opt.t)-obj.v(pipe))/obj.opt.t,pipe);
                end
            end
        end

        function obj=clean_buffer(obj)
            for pipe=obj.index.pipe.'
                obj.m_buffer{pipe}=[];
                obj.T_buffer{pipe}=[];
            end
        end

        function obj=add_m_buffer(obj,m,pipe)
            obj.m_buffer{pipe}=[obj.m_buffer{pipe} m];
        end

        function obj=add_T_buffer(obj,T,pipe)
            obj.T_buffer{pipe}=[obj.T_buffer{pipe} T];
        end

        function obj=update_T(obj, T)
            for pipe=obj.index.pipe.'
                if obj.direction(pipe)==1
                    obj.T{pipe}(1)=T(pipe);
                else
                    obj.T{pipe}(end)=T(pipe);
                end
            end
            obj.T_updated=true;
        end

        function obj=change_direction(obj, pipes)
            for pipe=reshape(pipes,1,[])
                obj.direction(pipe)=-obj.direction(pipe);
            end
        end
        
        function Tout=generate_Tout(obj)
            % generate outlet temperature of pipes
            Tout=zeros(obj.num,1);
            for pipe=obj.index.pipe.'
                [m_outlet,~]=get_outlet(obj,pipe);
                [m_inlet,~]=get_inlet(obj,pipe);
                Tout(pipe)=sum(obj.m_buffer{pipe}.*obj.T_buffer{pipe})/m_inlet;
                Tout(pipe)=obj.parameter.Ta+(Tout(pipe)-obj.parameter.Ta)*...
                    exp(-obj.parameter.lambda(pipe)/obj.parameter.area(pipe)/obj.parameter.density/obj.parameter.Cp*(length(obj.water_mass{pipe})-1+1/2+sum(obj.m_buffer{pipe}(2:end-1))/m_outlet)*obj.opt.t);
            end
        end

    end
end