function wx_up = upscale_update(wwx, image_level)


if image_level == 0
    wx_up = wwx;
else
    for ii = image_level:-1:1
        size_up = size(wwx)*2;
        
        [sx, sy, sz] = size(wwx);
        
        wx_up = single(zeros(size_up));
        
        wx_up(1:2:2*sx,1:2:2*sy,1:2:2*sz) = 2*wwx;
        
        
        wx_up(2:2:2*sx-1,1:2:2*sy,1:2:2*sz) = wwx(1:end-1,:,:) + wwx(2:end,:,:);
        
        wx_up(2*sx,1:2:2*sy,1:2:2*sz) = 2*wwx(end,:,:);
        
        
        
        wx_up(1:2:2*sx,2:2:2*sy-1,1:2:2*sz) = wwx(:,1:end-1,:) + wwx(:,2:end,:);
        
        wx_up(1:2:2*sx,2*sy,1:2:2*sz) = 2*wwx(:,end,:);
        
        
        
        wx_up(1:2:2*sx,1:2:2*sy,2:2:2*sz-1) = wwx(:,:,1:end-1) + wwx(:,:,2:end);
        
        wx_up(1:2:2*sx,1:2:2*sy,2*sz) = 2*wwx(:,:,end);
        
        
        
        wx_up(2:2:2*sx-1,2:2:2*sy,1:2:2*sz) = 0.5*(wx_up(1:2:2*sx-3,2:2:2*sy,1:2:2*sz)...
            + wx_up(3:2:2*sx,2:2:2*sy,1:2:2*sz));
        
        wx_up(2*sx,2:2:2*sy,1:2:2*sz) =  wx_up(2*sx-1,2:2:2*sy,1:2:2*sz);
        
        
        
        wx_up(1:2:2*sx,2:2:2*sy-1,2:2:2*sz) = 0.5*(wx_up(1:2:2*sx,1:2:2*sy-3,2:2:2*sz)...
            + wx_up(1:2:2*sx,3:2:2*sy,2:2:2*sz));
        
        wx_up(1:2:2*sx,2*sy,2:2:2*sz) =  wx_up(1:2:2*sx,2*sy-1,2:2:2*sz);
        
        
        
        wx_up(2:2:2*sx-1,1:2:2*sy,2:2:2*sz) = 0.5*(wx_up(1:2:2*sx-3,1:2:2*sy,2:2:2*sz)...
            + wx_up(3:2:2*sx,1:2:2*sy,2:2:2*sz));
        
        wx_up(2*sx,1:2:2*sy,2:2:2*sz) =  wx_up(2*sx-1,1:2:2*sy,2:2:2*sz);
        
        
        
        wx_up(2:2:2*sx-1,2:2:2*sy,2:2:2*sz) = 0.5*(wx_up(1:2:2*sx-3,2:2:2*sy,2:2:2*sz)...
            + wx_up(3:2:2*sx,2:2:2*sy,2:2:2*sz));
        
        
        wx_up(2*sx,2:2:2*sy,2:2:2*sz) = wx_up(2*sx-1,2:2:2*sy,2:2:2*sz);
        wwx = wx_up;
    end
end



return;