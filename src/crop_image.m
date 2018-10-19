% img_in:      2D image array input
% coordinates: 4 points in the img_in to be the edge points of the img_out;
% starting at the upper left vorner, the coordinates iterate clockwise
% through the image
% resolution:  the desired output resolution of the image [x,y]

function [img_out] = crop_image(img_in, coordinates_x, coordinates_y, resolution) 

debug_header = false;
debug_runtime = false;
debug_runtime_extra = false;
debug_polyshape = false;

% INCRESE RESOLUTION BY 1
    % it is increased by one such that there are coordinates for the pixels
    % to end properly; e.g. 0 to 5 in 10 steps means the first pixel goes
    % from 0.0 to 0.5, the second from 0.5 to 1.0 and the last one from 4.5
    % to 5.0; but keep in mind: 10 pixels mean that the last coordinate is
    % 4.5 but you need an endvalue for the calculation later on
    resolution = resolution + 1;

% CALCULATE X AND Y EDGE POINTS
    % get length
    len_12 = sqrt((coordinates_x(2) - coordinates_x(1))^2 + (coordinates_y(2) - coordinates_y(1))^2);
    len_34 = sqrt((coordinates_x(4) - coordinates_x(3))^2 + (coordinates_y(4) - coordinates_y(3))^2);
    len_23 = sqrt((coordinates_x(3) - coordinates_x(2))^2 + (coordinates_y(3) - coordinates_y(2))^2);
    len_41 = sqrt((coordinates_x(1) - coordinates_x(4))^2 + (coordinates_y(1) - coordinates_y(4))^2);
    
    if (debug_header)
        disp('LENGTH:')
        disp('- 12: ' + string(len_12));
        disp('- 23: ' + string(len_23));
        disp('- 34: ' + string(len_34));
        disp('- 41: ' + string(len_41));
        disp(' ');
    end
        
% CALCULATE X AND Y VALUES TO BE ADDED EACH SUBSTEP FOR EACH SIDE
    y_12 = (coordinates_y(2) - coordinates_y(1)) / (resolution(2) -1);
    x_12 = (coordinates_x(2) - coordinates_x(1)) / (resolution(1) -1);
    y_23 = (coordinates_y(3) - coordinates_y(2)) / (resolution(2) -1);
    x_23 = (coordinates_x(3) - coordinates_x(2)) / (resolution(1) -1);
    y_34 = (coordinates_y(4) - coordinates_y(3)) / (resolution(2) -1);
    x_34 = (coordinates_x(4) - coordinates_x(3)) / (resolution(1) -1);
    y_41 = (coordinates_y(1) - coordinates_y(4)) / (resolution(2) -1);
    x_41 = (coordinates_x(1) - coordinates_x(4)) / (resolution(1) -1);
    
    if (debug_header)
        disp('STEPS:');
        disp('- y_12: ' + string(y_12));
        disp('- x_12: ' + string(x_12));
        disp('- y_23: ' + string(y_23));
        disp('- x_23: ' + string(x_23));
        disp('- y_34: ' + string(y_34));
        disp('- x_34: ' + string(x_34));
        disp('- y_41: ' + string(y_41));
        disp('- x_41: ' + string(x_41));
        disp(' ');
    end
    
% CALCULATE OUTER PIXEL VALUES
    vec_y_12 = get_pixel_array(coordinates_y(1), y_12, resolution(1)); %resolution(1) because 12 represents the x axis of the image
    vec_x_12 = get_pixel_array(coordinates_x(1), x_12, resolution(1));
    vec_y_23 = get_pixel_array(coordinates_y(2), y_23, resolution(2));
    vec_x_23 = get_pixel_array(coordinates_x(2), x_23, resolution(2));
    vec_y_34 = get_pixel_array(coordinates_y(3), y_34, resolution(1));
    vec_x_34 = get_pixel_array(coordinates_x(3), x_34, resolution(1));
    vec_y_41 = get_pixel_array(coordinates_y(4), y_41, resolution(2));
    vec_x_41 = get_pixel_array(coordinates_x(4), x_41, resolution(2));
    
    if (debug_header)
        disp('SIDE PIXEL VECTORS:');
        disp('- vec_y_12: ' + strjoin(string(round(vec_y_12, 2)), {', '}));
        disp('- vec_x_12: ' + strjoin(string(round(vec_x_12, 2)), {', '}));
        disp('- vec_y_23: ' + strjoin(string(round(vec_y_23, 2)), {', '}));
        disp('- vec_x_23: ' + strjoin(string(round(vec_x_23, 2)), {', '}));
        disp('- vec_y_34: ' + strjoin(string(round(vec_y_34, 2)), {', '}));
        disp('- vec_x_34: ' + strjoin(string(round(vec_x_34, 2)), {', '}));
        disp('- vec_y_41: ' + strjoin(string(round(vec_y_41, 2)), {', '}));
        disp('- vec_x_41: ' + strjoin(string(round(vec_x_41, 2)), {', '}));
        disp(' ');
    end
    
% SWITCH VEC41 TO VEC14
    vec_length = resolution(2); %x and y are identical lengthwise

    vec_x_14 = zeros(vec_length, 1);
    vec_y_14 = zeros(vec_length, 1);
    for i = 1:1:vec_length
        vec_x_14(i) = vec_x_41(vec_length -i +1);
        vec_y_14(i) = vec_y_41(vec_length -i +1);
    end
    
    
% CALCULATE NEW PIXEL ARRAY
    pixels_x = zeros(resolution(2), resolution(1));
    pixels_y = zeros(resolution(2), resolution(1));
    
    for i = 1:1:resolution(2)
        A = get_coord_row(vec_x_14(i), vec_y_14(i), vec_x_23(i), vec_y_23(i), resolution(1));
        %disp(coordinates_y(4));
        %disp(vec_y_41(i));
        pixels_x(i,:) = A(1,:);
        pixels_y(i,:) = A(2,:);
    end
    
    if (debug_header)
        disp('COORDINATES X:');
        disp(pixels_x);
        disp(' ');
        disp('COORDINATES Y:');
        disp(pixels_y);
        disp(' ');
    end
    
    %output image
    img_out = zeros(resolution(2)-1, resolution(1)-1);
    
    if (debug_polyshape)
        fig1 = figure();
        figure(fig1); hold on;
        fig2 = figure();
        figure(fig2); hold on;
    end
    
% ITERATE THROUGH NEW PIXEL ARRAY
    for i=1:1:resolution(2)-1
        for k=1:1:resolution(1)-1
            % get vector of affected piexels
            min_x = floor(round(min([pixels_x(i,k), pixels_x(i,k+1), pixels_x(i+1,k), pixels_x(i+1,k+1)]),3));
            max_x = ceil(round(max([pixels_x(i,k), pixels_x(i,k+1), pixels_x(i+1,k), pixels_x(i+1,k+1)]),3) -1);
            min_y = floor(round(min([pixels_y(i,k), pixels_y(i,k+1), pixels_y(i+1,k), pixels_y(i+1,k+1)]),3));
            max_y = ceil(round(max([pixels_y(i,k), pixels_y(i,k+1), pixels_y(i+1,k), pixels_y(i+1,k+1)]),3) -1);
            
            if (debug_runtime_extra)
                disp('min_x: ' + string(min_x));
                disp('min_y:' + string(min_y));
                disp('max_x: ' + string(max_x));
                disp('max_y:' + string(max_y));
            end
            
            % new_pixel polyshape
            poly_new_pixel = polyshape(...
                [pixels_x(i,k), pixels_x(i,k+1), pixels_x(i+1,k+1), pixels_x(i+1,k)], ...
                [pixels_y(i,k), pixels_y(i,k+1), pixels_y(i+1,k+1), pixels_y(i+1,k)] ...
            ); %[x], [y]
            area_new_pixel = polyarea(poly_new_pixel.Vertices(:,1), poly_new_pixel.Vertices(:,2));
            
            if (debug_polyshape)
                figure(fig1);
                plot(poly_new_pixel);
            end
            
            if (debug_runtime)
                disp('new pixel area:' + string(area_new_pixel));
            end
            
            current_color_value = double(0.0);
            
            % calculate the overlap of the affected pixels
            for m=min_x:1:max_x
                for n=min_y:1:max_y
                    poly_old_pixel = polyshape([m,m,m+1,m+1], [n,n+1,n+1,n]);
                    poly_overlap = intersect(poly_new_pixel, poly_old_pixel);
                    area_overlap = polyarea(poly_overlap.Vertices(:,1), poly_overlap.Vertices(:,2));
                    
                    current_color_value = current_color_value + area_overlap * double(img_in(n+1,m+1));
                    
                    if (debug_polyshape)
                        figure(fig2);
                        %plot(poly_old_pixel);
                    end
                    
                    if (debug_runtime_extra)
                        disp('m: ' + string(m));
                        disp('n: ' + string(n));
                    end
                    
                    if (debug_runtime)
                        disp('overlap area: ' + string(area_overlap) + ', with colorvalue: ' + string(img_in(n+1,m+1)) + ' --> current color/area: ' + string(current_color_value));
                    end
               
                end
            end
            
            img_out(i,k) = current_color_value / area_new_pixel;
            
            if (debug_runtime)
                disp('new color value: ' + string(img_out(i,k)));
            end
            
            
            % iterate from mix_x to max_x and min_y to max_y
            % create 1x1 orig_pixel polyshape and overlap this with
            % new_pixel polyshape
            % overlap_area/new_pixel_area = percentage
            % colorvalue * percentage
            % add them up
            % -> new pixel           
            
            %poly1 = polyshape([0 0.5 1 0.5],[1 0 1 2]);
            %poly2 = polyshape([0.75 1.25 1.25 0.75],[0.25 0.25 0.75 0.75]);
            %plot(poly1)
            %hold on
            %plot(poly2)
            %disp(out.Vertices)
            %disp(polyarea(out.Vertices(:,1), out.Vertices(:,2)))
        end
        disp('Progress: ' + string(round(i/(resolution(2)-1)*100,2)) + '%');
    end
    
    
end

function [pixel_array] = get_pixel_array(startvalue, stepvalue, stepnumber)

    pixel_array = zeros(stepnumber, 1);
    pixel_array(1) = startvalue;
    for i = 2:1:stepnumber
        pixel_array(i) = pixel_array(i-1) + stepvalue;
    end

end

function [coord_row] = get_coord_row(start_x, start_y, end_x, end_y, stepnumber)

    %disp('start_x: ' + string(start_x) + ', end_x: ' + string(end_x));
    %disp('start_y: ' + string(start_y) + ', end_y: ' + string(end_y));

    coord_row_x = zeros(1, stepnumber);
    coord_row_y = zeros(1, stepnumber);
    step_x = (end_x - start_x) / (stepnumber -1);
    step_y = (end_y - start_y) / (stepnumber -1);
    
    coord_row_x(1) = start_x;
    coord_row_y(1) = start_y;
    
    for i = 2:1:stepnumber 
        coord_row_x(i) = coord_row_x(i-1) + step_x;
        coord_row_y(i) = coord_row_y(i-1) + step_y;
    end
    
    coord_row = [coord_row_x; coord_row_y];

end