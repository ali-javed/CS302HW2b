function comparederivs(allx,f,truedf)
%testing github
% COMPAREDERIVS: % Compares true, 1st-order, and 2nd-order slope approximations
% function comparederivs(allx,f,truedf)
%
% INPUTS:
%   allx: is a scalar or vector of what domain values to look at
%   f: is a function handle of the function to differentiate
%   truedf: is a function handle for the analytic derivative function of f
%
% OUTPUTS:
%   This program produces one plot for each value in the vector allx
%   The plots show the error in the approximations as a function of step size
%
% sample call: comparederivs([0 5],@exp, @exp)

% Maggie Eppstein, 2/10/08; documentation improved 2/15/11

% THIS IS A GOOD EXAMPLE OF A WELL-DOCUMENTED FILE;  NOTE THE FOLLOWING:
%   a) contents and consistent organization of function headers; first
%   comment line for LOOKFOR command, 1st contiguous comment block for HELP
%   command; always define inputs/outputs, including size constraints or
%   other pre-/post-conditions;
%   b) in-line comments should always be at one level of abstraction higher
%   than the code itself;
%   c) use of full-line UPPER-CASE in-line comments to give a high-level
%   description of what each logically-related code block does; you can read through these alone to
%   get a good understanding of what the code does, without even looking at
%   the code;
%   d) additional lower-case comments at the ends of potentially confusing
%   lines for clarification;


% FOR EACH X-VALUE SPECIFIED, PLOT THE APPROXIMATION ERRORS AS A FUNCTION OF STEPSIZE
h=logspace(-20,-1,20); %logarithmically-spaced step sizes to try
for xi=1:length(allx) %each x gets its own plot
    x=allx(xi);
    
    % COMPUTE TRUE DERIVATIVE AND ITS APPROXIMATIONS
    df=truedf(x); % compute true derivative at x
    df1=backdiff(x,f,h); %approximate with 1st order backwards difference
    df2=centraldiff(x,f,h); %approximate with 2nd order central difference
    
    % COMPUTE APPROXIMATION ERRORS BY COMPARING TO TRUE DERIVATIVES
    err1=abs(df-df1);
    err2=abs(df-df2);

    % PLOT THE APPROXIMATION ERRORS AS A FUNCTION OF STEPSIZE
    figure
    loglog(h,err1,'ro',h,err2,'b*');
    
    % LINEAR REGRESSION OF LOG-LOG RELATIONSHIPS (ONLY IN REGION GOVERNED BY TRUNCATION ERROR)
    coef1=polyfit(log(h(end-4:end)),log(err1(end-4:end)),1);
    coef2=polyfit(log(h(end-4:end)),log(err2(end-4:end)),1);
    
    hold on % add the best-fit lines to the plot
    loglog(h(end-4:end),exp(coef1(2))*h(end-4:end).^coef1(1),'r-');
    loglog(h(end-4:end),exp(coef2(2))*h(end-4:end).^coef2(1),'b-');
    
    % LABEL THE PLOT
    set(gca,'fontsize',14) % be kind to the instructor's aging eyes!
    xlabel('h')
    ylabel('error')
    legend('back diff','central diff',...
        ['slope = ',num2str(coef1(1),4)],['slope = ',num2str(coef2(1),4)],...
        'Location','BestOutside');
    title(['at x = ',num2str(x)])
    
    %ALLOW USER TO VIEW EACH PLOT BEFORE MOVING ON TO THE NEXT
    if xi<length(allx)
        disp('Hit any key to continue...')
        pause 
    end
end
figure(gcf)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: the following functions are placed here for convenience for this demo code, 
% but cannot be called from outside this file; in general, these should be
% in their own files (e.g. in a directory for your personal library "toolbox" that 
% you add to the Matlab path to access your own handy utility functions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function df1=backdiff(x,f,h) 
% BACKDIFF: 1st order backwards difference approximation to first derivative
% function df1=backdiff(x,f,h) 
%
% INPUTS:
% x: location(s) of where in domain to approximat the derivative
% f: handle of function to approximate derivative of
% h: stepsize(s) to use in approximation
% SIZE CONSTRAINTS: at least one of x or h must be a scalar, but the other
% can be of any dimension (scalar, vector, matrix)
%
% OUTPUTS:
% df1: 1st order approximation to first deriv (slope) of f at x 
%     (same size as largest of x or h)
% 
% SAMPLE CALLS:
%   df1=backdiff(0,@sin,[1e-3 1e-2 1e-1]) % can call with vector of stepsizes
%   df1=backdiff(0:.5:3,@sin,1e-3) % or can call with vector of domain values

% AUTHOR: Maggie Eppstein, 2/15/2011

df1=(f(x)-f(x-h))./h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df2=centraldiff(x,f,h)
% CENTRALDIFF: 2nd order central difference approximation to first derivative
%
% INPUTS:
% x: location(s) of where in domain to approximat the derivative
% f: handle of function to approximate derivative of
% h: stepsize(s) to use in approximation
% SIZE CONSTRAINTS: at least one of x or h must be a scalar, but the other
% can be of any dimension (scalar, vector, matrix)
%
% OUTPUTS:
% df2: 2nd order approximation to first deriv (slope) of f at x 
%     (same size as largest of x or h)
%
% SAMPLE CALLS:
%   df2=backdiff(0,@sin,[1e-3 1e-2 1e-1]) % can call with vector of stepsizes
%   df2=backdiff(0:.5:3,@sin,1e-3) % or can call with vector of domain values

% AUTHOR: Maggie Eppstein, 2/15/2011

df2=(f(x+h)-f(x-h))./(2*h);