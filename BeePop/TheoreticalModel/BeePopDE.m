classdef BeePopDE < handle
    properties 
        tend;
        dt;
        ic;
        t;
        beta  = 1;
        gamma = 1;
        mu    = 1; 
        delta = @(r) 1;
        dDdr  = @(r) 0;
        result;
        lam;
        V;
        Frames;
        Population;
    end
    methods
        function bp = BeePopDE(tend, dt, ic, pop)
            bp.tend  = tend;
            bp.dt    = dt;
            bp.t     = 0:dt:bp.tend;
            bp.ic    = ic;
            bp.Population = pop;
        end
        function Simulate(obj, isPeriodic)
            if isPeriodic
                [~,y] = ode45(@obj.flowPeriodic,obj.t,obj.ic);
            else
                [~,y] = ode45(@obj.flow,obj.t,obj.ic);
            end
        obj.result = y;
        end
        function dy = Update(obj,y)
            dy = y + obj.flow(obj.t,y)*obj.dt;
        end
        function dy = flow(obj,~, y)
            
            dy    =  zeros(size(y));   
            dy(1) = -obj.beta*y(1)*y(2);
            dy(2) =  obj.beta*y(1)*y(2) - obj.gamma*y(2);
            dy(3) =  obj.gamma*y(2);
        end      
        function dy = flowPeriodic(obj,~, y)
            dy    =  zeros(size(y));   
            %dy(1) = -obj.beta*y(1)*y(2) + obj.mu*(1-y(1));
            %dy(2) =  obj.beta*y(1)*y(2) - (obj.gamma + obj.mu)*y(2);
            %dy(3) =  obj.gamma*y(2) - (obj.mu)*y(3);
            dy(1) = -obj.beta*y(1)*y(2) + obj.mu*y(3);
            dy(2) =  obj.beta*y(1)*y(2) - obj.gamma*y(2);
            dy(3) =  obj.gamma*y(2) - obj.mu*y(3);
        end       
        function J = Jacobian(obj, t0)
            id = find(t0-obj.dt/2 <= obj.t & obj.t < t0+obj.dt/2); 
            S  = obj.result(id,1);
            I  = obj.result(id,2);
            R  = obj.result(id,3);
            a  = obj.beta;
            b  = obj.gamma;
            
            J = [-a*I*obj.delta(R)*S^(obj.delta(R)-1), ... 
                 -a*S^obj.delta(R),...
                 -a*I*log(S)*S^obj.delta(R)*obj.dDdr(R);...
                 
                 a*I*obj.delta(R)*S^(obj.delta(R)-1),...
                 a*S^obj.delta(R) - b,...
                 a*I*S^obj.delta(R)*log(S)*obj.dDdr(R);...
                 
                 0,...
                 b,...
                 0];
        end 
        function getEigen(obj)
            obj.lam = zeros(3,length(obj.t));
            obj.V = zeros(3,3,length(obj.t));
            for ii = 1:length(obj.t)
                J = obj.Jacobian(obj.t(ii));
                [v,l] = eig(J);
                obj.lam(:,ii) = diag(l);
                obj.V(:,:,ii) = v;
                obj.lam(abs(obj.lam)<1e-5)=0;
            end
        end
        function plot(obj,leg,varargin)
            S = obj.result(:,1)*obj.Population;
            I = obj.result(:,2)*obj.Population;
            R = obj.result(:,3)*obj.Population;
            %figure
            plot(obj.t, S, 'b--',obj.t, I, 'r--',obj.t, R, 'g--','LineWidth',3,varargin{:})
            if strcmpi(leg,'LegendOn')
                l=legend('Susceptible','Infected','Recovered');
                set(l,'FontSize',18)
            end
            ti = title(sprintf('\\beta = %3.2e, \\gamma = %3.2e',obj.beta, obj.gamma));
            set(ti, 'FontSize',20)
            xlabel('Time','FontSize',18)
            ylabel('S,I,R','FontSize',18)
        end    
        function timePlotEigen(obj)
            lam1 = obj.lam(1,:);
            lam2 = obj.lam(2,:);
            lam3 = obj.lam(3,:);
            
            figure
            subplot(3,1,1)
            plot(obj.t,real(lam1),'k','LineWidth',3)
            xlabel('Time','FontSize',18)
            ylabel('\lambda_1','FontSize',18)
            ylim([-1,1])
            hold on
            subplot(3,1,2)
            plot(obj.t,real(lam2),'k','LineWidth',3)
            xlabel('Time','FontSize',18)
            ylabel('\lambda_2','FontSize',18)
            hold on
            subplot(3,1,3)
            plot(obj.t,real(lam3),'k','LineWidth',3)
            xlabel('Time','FontSize',18)
            ylabel('\lambda_3','FontSize',18)
            hold on
        end
        function argonPlotEigen(obj)
            lam1 = obj.lam(1,:);
            lam2 = obj.lam(2,:);
            lam3 = obj.lam(3,:);
            
            figure
            subplot(2,1,1)
            
            plot(real(lam1),imag(lam1),'k.')
            title('\lambda_1','FontSize',20)
            xlabel('Real Part','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
           
            subplot(2,1,2)
            plot(real(lam2),imag(lam2),'k.',...
                 real(lam3),imag(lam3),'b.')
            title('\lambda_2 and \lambda_3','FontSize',20)
            xlabel('Real Part','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
            l = legend('\lambda_2','\lambda_3');
            set(l,'FontSize',18)
        end
        function plotEigenVector(obj,t)
            v = obj.V(:,:,t);
            figure
            for ii = 1:length(v)
                quiver3(0,0,0,v(ii,1),v(ii,2),v(ii,3))
                hold on;
            end
            hold off
        end
        function h = singlePlot(obj,idx, fignum, varargin)
            
            if idx == 'All'
                idx = length(obj.t);
            end
            
            S = obj.result(1:idx,1);
            I = obj.result(1:idx,2);
            R = obj.result(1:idx,3);
            
            lam1 = obj.lam(1,1:idx);
            lam2 = obj.lam(2,1:idx);
            lam3 = obj.lam(3,1:idx);
            
            time = obj.t(1:idx);
            
            h = figure(fignum);
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            subplot(3,2,[1,2])
            
            plot(time, S, time, I, time, R,'LineWidth',3,varargin{:})
            grid on
            axis([0 obj.tend 0 1]);
            l=legend('Susceptible','Infected','Recovered','Location','northoutside','Orientation','horizontal');
            set(l,'FontSize',18)
            
            %ti = title(sprintf('\\beta = %3.2e, \\gamma = %3.2e',obj.beta, obj.gamma));
            %set(ti, 'FontSize',20)
            xlabel('Time','FontSize',18)
            ylabel('S,I,R','FontSize',18)
            
            subplot(3,2,[3,5])
            plot(real(lam2),imag(lam2),'k.',...
                 real(lam3),imag(lam3),'r.')
            grid on

            axis([min(min(real(obj.lam))), max(max(real(obj.lam))),...
                 min(min(imag(obj.lam))), max(max(imag(obj.lam)))])
            
            hold on
            line('XData', [min(min(real(obj.lam))), max(max(real(obj.lam)))],...
                 'YData', [0, 0]);
            hold on
            line('XData', [0, 0],...
                 'YData', [min(min(imag(obj.lam))), max(max(imag(obj.lam)))]);
             
            xlabel('Real Part','FontSize',18)
            ylabel('Imaginary Part','FontSize',18)
            l = legend('\lambda_2','\lambda_3','Location','northeast');
            set(l,'FontSize',18)
             
            
            subplot(3,2,4)
            plot(time,real(lam2),'k','LineWidth',3)
            grid on
            hold on
            line('XData', [0,obj.tend], 'YData',[0,0]);
            axis([0, obj.tend, min(min(real(obj.lam))),  max(max(real(obj.lam)))]);
            xlabel('Time','FontSize',18)
            ylabel('\lambda_2','FontSize',18)
            subplot(3,2,6)
            plot(time,real(lam3),'k','LineWidth',3)
            grid on
            hold on
            line('XData', [0,obj.tend], 'YData',[0,0]);
            axis([0, obj.tend, min(min(real(obj.lam))),  max(max(real(obj.lam)))]);
            xlabel('Time','FontSize',18)
            ylabel('\lambda_3','FontSize',18)
            
        end
        function makeMovie(obj,fignum,varargin)
            frames(length(obj.t)) = struct('cdata',[],'colormap',[]);
            for ii = 1:length(obj.t)
                h = obj.singlePlot(ii,fignum);
                frames(ii) = getframe(h);
                close(h)
            end   
            obj.Frames = frames;
        end
    end
end
