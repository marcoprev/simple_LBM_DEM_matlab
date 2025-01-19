close all
clearvars
% Domain parameters
Nx = 400;
Ny = 100;
mean_rho = 100;  % average density (non dim)
% Bhatnagar-Gross-Krook (BGK) model: Omega(f)=-1/tau*(f-feq)
tau = 0.6;   % relaxation time. viscosity = speed of sound ^2 *(tau-deltaT/2)
% speed of sound depends on the discretisation adopted (e.g. 1/sqrt(3)dxdt)

dt = 2e-1;   % s
nsteps = 1000;   

% Lattice weights (D2Q9)
NL = 9;
% velocity sets
vec_cx = [0, 0, 1, 1, 1, 0, -1, -1, -1];% * dx/dt
vec_cy = [0, 1, 1, 0, -1, -1, -1, 0, 1];
vec_weights = [4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36];

% sum of each discrete population for each lattice direction
%sum(vec_c(i)*mat_F(i))=rho*mat_u

    % Initial Conditions
    mat_F = ones(Ny, Nx, NL); % Constant distribution

    [mat_X, mat_Y] = meshgrid(1:Nx, 1:Ny);
    % set the velocity by modifying the horizontal component (4th direction) on the F field 
    mat_F(:, :, 4) = mat_F(:, :, 4) + 2;

    % scale mat_F so the total rho is equal to the target one
    mat_rho = sum(mat_F, 3);
    for i = 1:NL
        mat_F(:, :, i) = mat_F(:, :, i) .* (mean_rho ./ mat_rho);
    end

% outputVideo = VideoWriter('out_lbm.mp4','MPEG-4');
% framerate = 10;
% outputVideo.FrameRate = framerate;
% open(outputVideo)
% 


    % Create a particle as a cylinder boundary in 2D
    cylR = 10;
    cylMass = cylR^2*pi*mean_rho*1.5;

    % initial pos and velocity
    cylPos = [cylR+10,Ny/2];
    x0 = cylPos(1);
    cylVel =[0,0];


   for step_idx = 1:nsteps

        % apply periodic conditions in the horizontal direction
        for i = 1:NL
            mat_F(:, :, i) = circshift(mat_F(:, :, i), [vec_cy(i), vec_cx(i)]);
        end

        % update ball position
        cylPos=cylPos+cylVel*dt;

        % apply the particle as boundary to the fluid (staircase approach)
        mat_cylinder = sqrt((mat_X - cylPos(1)).^2 + (mat_Y - cylPos(2)).^2) < cylR;


        % use it to set reflective boundaries, i.e. invert the direction
        ball_idx = find(mat_cylinder); % identify cells part of the ball
        mat_F2D = reshape(mat_F, [], NL); % big matrix of all directions in one
        mat_F_ball = mat_F2D(ball_idx, :); % get the F value in this big matrix at the ball pos
        mat_F_ball_reflected = mat_F_ball(:, [1, 6, 7, 8, 9, 2, 3, 4, 5]); % invert them


        % use this difference to get the total change in F
        mat_deltaF = mat_F_ball_reflected - mat_F_ball;
        vec_fx = sum(mat_deltaF .* vec_cx, 2); % mult by its normals
        vec_fy = sum(mat_deltaF .* vec_cy, 2);
        forceX = sum(vec_fx);
        forceY = sum(vec_fy);
       
        cylAcc = -[forceX,forceY]/cylMass;
        cylVel = cylVel + cylAcc*dt;
        
        % re-apply them to the current field
        for i = 1:NL
            tmp_F = mat_F(:, :, i);
            tmp_F(mat_cylinder) = mat_F_ball_reflected(:, i);
            mat_F(:, :, i) = tmp_F;
        end

        % retreive fluid variables (density and velocities)
        mat_rho = sum(mat_F, 3);
        mat_ux = sum(mat_F .* reshape(vec_cx, 1, 1, []), 3) ./ mat_rho;
        mat_uy = sum(mat_F .* reshape(vec_cy, 1, 1, []), 3) ./ mat_rho;

        % use some operator to apply particle collisions
        mat_Feq = zeros(size(mat_F));
        for i = 1:NL
            cx = vec_cx(i);
            cy = vec_cy(i);
            w = vec_weights(i);
            % the collision operator
            % feq = w * rho * ( 1 + 3*c*u + 9/2*(c*u)^2 - 3/2 u^2;
            mat_Feq(:, :, i) = mat_rho .* w .* (1 + 3 * (cx * mat_ux + cy * mat_uy) + ...
                9/2 * (cx * mat_ux + cy * mat_uy).^2 - 3/2 * (mat_ux.^2 + mat_uy.^2));
        end
        % update the F field depending on delta F and tau
        mat_F = mat_F - (1 / tau) * (mat_F - mat_Feq);

            % set 0 velocities at the ball location
            mat_ux(mat_cylinder) = 0;
            mat_uy(mat_cylinder) = 0;
           
            clf
            imagesc(sqrt(mat_ux.^2+mat_uy.^2), 'AlphaData', ~mat_cylinder);
            viscircles(cylPos,cylR,'Color','k');


            hold on

            %subsample the matrices to plot vectors
        myScale = 0.3;
          quiver(imresize(mat_X,myScale),imresize(mat_Y,myScale),imresize(mat_ux,myScale),imresize(mat_uy,myScale),'Color','k');
          colorbar
       
          title(['Time: ' num2str(dt*step_idx) ', ball displacement:' num2str(cylPos(1)-x0)])
            axis equal;
            drawnow;
% frame = getframe(gcf);
% im = frame2im(frame);
%    writeVideo(outputVideo,im)
% 

        
    end

%close(outputVideo)
