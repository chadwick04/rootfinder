function root_finder_ui
% Simple UI for six root-finding methods with auto root detection in Graphical
% Added dynamic instructions per method

%% ---------- FIGURE & BASIC CONTROLS ------------------------------------
fig = uifigure('Name','Numerical Methods Root Finder',...
               'Position',[100 100 820 680]);

% — Equation —
uilabel(fig,'Text','Equation f(x):','Position',[20 615 100 22]);
eqField = uieditfield(fig,'text','Position',[125 615 350 22], ...
        'Value','x.^3 - 6*x.^2 + 11*x - 6');       % default demo

% — Method —
uilabel(fig,'Text','Method:','Position',[20 575 100 22]);
methodDrop = uidropdown(fig,'Items',...
       {'Graphical','Incremental Search','Bisection',...
        'False Position','Newton-Raphson','Secant'},...
       'Position',[125 575 200 22],...
       'ValueChangedFcn',@(dd,~) updateInputs(dd.Value));

% — Parameter fields (2 generic numeric boxes) —
param1Label = uilabel(fig,'Position',[20 535 140 22]);
param1Field = uieditfield(fig,'numeric','Position',[165 535 120 22]);
param2Label = uilabel(fig,'Position',[295 535 140 22]);
param2Field = uieditfield(fig,'numeric','Position',[440 535 120 22]);

% — Instructions TextArea (multiline) —
instrArea = uitextarea(fig,...
    'Position',[20 430 700 80],...
    'Editable',false,...
    'FontName','Consolas',...
    'FontSize',11);

% — Solve button —
uibutton(fig,'Text','Solve','Position',[600 525 100 30],...
         'ButtonPushedFcn',@(btn,~) solve());

% — Axes & data table —
ax  = uiaxes(fig,'Position',[50 170 720 250]);
tbl = uitable(fig,'Position',[50 20 720 130]);

%% ---------- INITIAL UI STATE -------------------------------------------
updateInputs(methodDrop.Value)                 % set labels & instructions for default

%% ---------- NESTED HELPER FUNCTIONS ------------------------------------
    function updateInputs(method)
        % Reset fields
        param1Field.Value = 0;  param2Field.Value = 0;
        param2Field.Visible = 'on';  param2Label.Visible = 'on';
        switch method
            case {'Bisection','False Position'}
                param1Label.Text = 'Lower bound  a:';         % a
                param2Label.Text = 'Upper bound  b:';         % b
                instrArea.Value = {
                    'Bisection / False Position Method';
                    '---------------------------------';
                    'Purpose: Find a root by narrowing interval [a,b]';
                    'Usage: Enter a and b where f(a)*f(b)<0 (sign change)';
                    'Bisection: root is midpoint of interval each iteration';
                    'False Position: root estimate uses linear interpolation';
                    'Parameters: Lower bound a, Upper bound b';
                    'Output: Approximate root near [a,b]';
                    'Notes: Converges reliably if sign change exists';
                    };
            case 'Incremental Search'
                param1Label.Text = 'Start x:';                % x_start
                param2Label.Text = 'dx (step):';              % delta
                instrArea.Value = {
                    'Incremental Search Method';
                    '-------------------------';
                    'Purpose: Find intervals containing roots by stepping';
                    'Usage: Enter starting x and step size dx';
                    'The method checks sign changes at increments of dx';
                    'Parameters: Start x, step dx (must be >0)';
                    'Output: Approximate roots where sign changes occur';
                    'Notes: Step size affects root detection accuracy';
                    };
            case 'Newton-Raphson'
                param1Label.Text = 'Initial guess  x0:';      % x0
                param2Label.Text = '';                        % hide 2nd
                param2Field.Visible = 'off'; param2Label.Visible = 'off';
                instrArea.Value = {
                    'Newton-Raphson Method';
                    '----------------------';
                    'Purpose: Iterative root finding using derivative';
                    'Usage: Enter initial guess x0';
                    'Parameters: Initial guess x0';
                    'Output: Approximate root from iteration';
                    'Notes: Requires differentiable function';
                    'May fail if derivative near zero or poor guess';
                    };
            case 'Secant'
                param1Label.Text = 'First guess  x0:';        % x0
                param2Label.Text = 'Second guess  x1:';       % x1
                instrArea.Value = {
                    'Secant Method';
                    '--------------';
                    'Purpose: Iterative root finding using two initial guesses';
                    'Usage: Enter two initial guesses x0 and x1';
                    'Parameters: First guess x0, second guess x1';
                    'Output: Approximate root from iteration';
                    'Notes: Does not require derivative explicitly';
                    };
            otherwise   % Graphical
                param1Label.Text = '';  param2Label.Text = '';
                param2Field.Visible = 'off'; param2Label.Visible = 'off';
                instrArea.Value = {
                    'Graphical Method';
                    '------------------';
                    'Purpose: Visualize function and estimate roots';
                    'Usage: No parameters needed';
                    'Method: Plots f(x) over range [-10,10]';
                    'Roots are detected by sign changes and marked';
                    'Output: Approximate roots shown as red circles';
                    'Notes: Useful for initial guess and multiple roots';
                    };
        end
    end

%-------------------------------------------------------------------------    
    function solve()
        cla(ax); tbl.Data = {};                % clear old
        f_str = eqField.Value;
        if isempty(f_str), uialert(fig,'Enter a function!','Input Error'); return; end
        try
            f = str2func(['@(x) ' f_str]);
            % Test function at 0 for early error catch
            f(0);
        catch
            uialert(fig,'Invalid function expression!','Input Error');
            return;
        end
        
        method = methodDrop.Value;
        switch method
        % ---------------- GRAPHICAL --------------------------------------
        case 'Graphical'
            % Fine grid for detecting roots automatically
            x = linspace(-10,10,10000);  
            y = f(x);
            plot(ax,x,y,'LineWidth',1.4); grid(ax,'on'); hold(ax,'on');
            plot(ax,x,zeros(size(x)),'k--');

            % Find sign changes => root intervals
            signChangeIdx = find(diff(sign(y))~=0);
            rootsFound = zeros(size(signChangeIdx));
            for k=1:length(signChangeIdx)
                idx = signChangeIdx(k);
                % Linear interpolation for root approx
                x1 = x(idx); x2 = x(idx+1);
                y1 = y(idx); y2 = y(idx+1);
                rootApprox = x1 - y1*(x2-x1)/(y2-y1);
                rootsFound(k) = rootApprox;
                plot(ax,rootApprox,0,'ro','MarkerSize',8,'LineWidth',1.5);
            end

            % Legend with roots numbered
            legendStr = ["f(x)","y=0"];
            for k=1:length(rootsFound)
                legendStr(end+1) = sprintf("Root%d: %.6g", k, rootsFound(k));
            end
            legend(ax,legendStr,'Location','best');

            xlabel(ax,'x'); ylabel(ax,'f(x)');
            title(ax,'Graphical Method (Automatic root detection)');

            tbl.ColumnName = {'Root Approximation'};
            tbl.Data = rootsFound';

        % --------------- INCREMENTAL SEARCH ------------------------------
        case 'Incremental Search'
            x0 = param1Field.Value;  dx = param2Field.Value;
            if dx<=0, uialert(fig,'dx must be >0','Input Error'); return; end
            xmax = x0 + 100*dx;             % scan 100 steps
            xs = x0:dx:xmax; ys = f(xs);
            plot(ax,xs,ys,'b-',xs,zeros(size(xs)),'k--'); hold(ax,'on');
            xlabel(ax,'x'); ylabel(ax,'f(x)'); grid(ax,'on');
            title(ax,'Incremental Search');
            roots = [];
            for i=1:numel(xs)-1
                if ys(i)*ys(i+1)<=0
                    r = (xs(i)+xs(i+1))/2;
                    roots(end+1) = r; %#ok<AGROW>
                    plot(ax,r,0,'ms','MarkerSize',10,'LineWidth',1.5);
                end
            end
            tbl.ColumnName = {'Root (mid-point)'};
            tbl.Data = roots';

        % --------------- BISECTION / FALSE POSITION ----------------------
        case {'Bisection','False Position'}
            a = param1Field.Value;  b = param2Field.Value;
            if a>=b, uialert(fig,'Require a<b','Input Error'); return; end
            fa = f(a); fb = f(b);
            if fa*fb>0
                uialert(fig,'f(a) & f(b) must have opposite sign','Error');
                return;
            end
            tol = 1e-6; maxIt = 50;
            iterData = [];
            for it=1:maxIt
                if strcmp(method,'Bisection')
                    c = (a+b)/2;
                else       % False Position
                    c = b - fb*(a-b)/(fa-fb);
                end
                fc = f(c);
                iterData(end+1,:) = [it, a, b, c, fc]; %#ok<AGROW>
                if abs(fc)<tol || abs(b-a)<tol, break; end
                if fa*fc<0, b=c; fb=fc; else, a=c; fa=fc; end
            end
            % Plot
            xx = linspace(a,b,400); yy = f(xx);
            plot(ax,xx,yy,'b-',xx,0*xx,'k--'); hold(ax,'on');
            plot(ax,c,0,'gd','MarkerSize',10,'LineWidth',1.5);
            xlabel(ax,'x'); ylabel(ax,'f(x)'); grid(ax,'on');
            title(ax,[method ' Method']);
            text(ax,c,0,sprintf('  %.6g',c),'Color','g','FontWeight','bold');
            tbl.ColumnName = {'Iter','a','b','c','f(c)'};
            tbl.Data = iterData;

        % --------------- NEWTON RAPHSON ----------------------------------
        case 'Newton-Raphson'
            x0 = param1Field.Value;
            syms x; df = matlabFunction(diff(str2sym(f_str)));
            tol=1e-6; maxIt=50;
            iterData=[];
            for it=1:maxIt
                fx = f(x0); dfx = df(x0);
                if dfx==0, break; end
                x1 = x0 - fx/dfx;
                iterData(end+1,:) = [it, x0, fx]; %#ok<AGROW>
                if abs(x1-x0)<tol, break; end
                x0 = x1;
            end
            root = x1;
            % Plot local view
            xp = linspace(root-2,root+2,400); yp=f(xp);
            plot(ax,xp,yp,'b-',xp,0*xp,'k--'); hold(ax,'on');
            plot(ax,root,0,'c^','MarkerSize',10,'LineWidth',1.5);
            xlabel(ax,'x'); ylabel(ax,'f(x)'); grid(ax,'on');
            title(ax,'Newton–Raphson');
            text(ax,root,0,sprintf('  %.6g',root),'Color','c','FontWeight','bold');
            tbl.ColumnName={'Iter','x','f(x)'};
            tbl.Data = iterData;

        % --------------- SECANT ------------------------------------------
        otherwise   % 'Secant'
            x0 = param1Field.Value;  x1 = param2Field.Value;
            tol=1e-6; maxIt=50; iterData=[];
            for it=1:maxIt
                f0=f(x0); f1=f(x1);
                if f1-f0==0, break; end
                x2 = x1 - f1*(x1-x0)/(f1-f0);
                iterData(end+1,:)=[it,x0,x1,x2,f(x2)]; %#ok<AGROW>
                if abs(x2-x1)<tol, break; end
                x0=x1; x1=x2;
            end
            root=x2;
            xp=linspace(root-2,root+2,400); yp=f(xp);
            plot(ax,xp,yp,'b-',xp,0*xp,'k--'); hold(ax,'on');
            plot(ax,root,0,'m*','MarkerSize',10,'LineWidth',1.5);
            xlabel(ax,'x'); ylabel(ax,'f(x)'); grid(ax,'on');
            title(ax,'Secant Method');
            text(ax,root,0,sprintf('  %.6g',root),'Color','m','FontWeight','bold');
            tbl.ColumnName={'Iter','x0','x1','x2','f(x2)'};
            tbl.Data=iterData;
        end
    end
end
