mergeInto(LibraryManager.library,{
    emcctime: function() {
	return Math.floor(Date.now());
    },
    glcontext:function(handle){
	GL.makeContextCurrent(handle);
    },
    glinit: function(width, height, depth, flags,no){
	var canvas;
	if (no==-1)
	    canvas=Module['canvas'];
	else
	    canvas=document.getElementById('gl3d_'+(no+1));

	// (0,0) means 'use fullscreen' in native; in Emscripten, use the current canvas size.
	if (width == 0 && height == 0) {
	    width = canvas.width;
	    height = canvas.height;
	}
	
	if (!SDL.addedResizeListener) {
	    SDL.addedResizeListener = true;
	    Browser.resizeListeners.push(function(w, h) {
		if (!SDL.settingVideoMode) {
		    SDL.receiveEvent({
			type: 'resize',
			w: w,
			h: h
		    });
		}
	    });
	}
	
	if (width !== canvas.width || height !== canvas.height) {
	    SDL.settingVideoMode = true; // SetVideoMode itself should not trigger resize events
	    Browser.setCanvasSize(width, height);
	    SDL.settingVideoMode = false;
	}
	
	// Free the old surface first if there is one
	if (SDL.screen) { SDL.freeSurface(SDL.screen); assert(!SDL.screen);	}
	
	if (SDL.GL) flags = flags | 0x04000000; // SDL_OPENGL - if we are using GL, then later calls to SetVideoMode may not mention GL, but we do need it. Once in GL mode, we never leave it.
	
      flags = flags || 0;
      var is_SDL_HWSURFACE = flags & 0x00000001;
      var is_SDL_HWPALETTE = flags & 0x00200000;
      var is_SDL_OPENGL = flags & 0x04000000;

      var surf = _malloc({{{ C_STRUCTS.SDL_Surface.__size__ }}});
      var pixelFormat = _malloc({{{ C_STRUCTS.SDL_PixelFormat.__size__ }}});
      //surface with SDL_HWPALETTE flag is 8bpp surface (1 byte)
      var bpp = is_SDL_HWPALETTE ? 1 : 4;
      var buffer = 0;

      // preemptively initialize this for software surfaces,
      // otherwise it will be lazily initialized inside of SDL_LockSurface
      if (!is_SDL_HWSURFACE && !is_SDL_OPENGL) {
        buffer = _malloc(width * height * 4);
      }

      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.flags, 'flags', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.format, 'pixelFormat', 'void*') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.w, 'width', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.h, 'height', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.pitch, 'width * bpp', 'i32') }}};  // assuming RGBA or indexed for now,
                                                                                        // since that is what ImageData gives us in browsers
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.pixels, 'buffer', 'void*') }}};

      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.clip_rect+C_STRUCTS.SDL_Rect.x, '0', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.clip_rect+C_STRUCTS.SDL_Rect.y, '0', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.clip_rect+C_STRUCTS.SDL_Rect.w, 'Module["canvas"].width', 'i32') }}};
      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.clip_rect+C_STRUCTS.SDL_Rect.h, 'Module["canvas"].height', 'i32') }}};

      {{{ makeSetValue('surf', C_STRUCTS.SDL_Surface.refcount, '1', 'i32') }}};

      {{{ makeSetValue('pixelFormat', C_STRUCTS.SDL_PixelFormat.format, cDefine('SDL_PIXELFORMAT_RGBA8888'), 'i32') }}};
      {{{ makeSetValue('pixelFormat', C_STRUCTS.SDL_PixelFormat.palette, '0', 'i32') }}};// TODO
      {{{ makeSetValue('pixelFormat', C_STRUCTS.SDL_PixelFormat.BitsPerPixel, 'bpp * 8', 'i8') }}};
      {{{ makeSetValue('pixelFormat', C_STRUCTS.SDL_PixelFormat.BytesPerPixel, 'bpp', 'i8') }}};


      // Decide if we want to use WebGL or not
      SDL.GL = SDL.GL || is_SDL_OPENGL;

      var webGLContextAttributes = {
        antialias: ((SDL.glAttributes[13 /*SDL_GL_MULTISAMPLEBUFFERS*/] != 0) && (SDL.glAttributes[14 /*SDL_GL_MULTISAMPLESAMPLES*/] > 1)),
        depth: (SDL.glAttributes[6 /*SDL_GL_DEPTH_SIZE*/] > 0),
        stencil: (SDL.glAttributes[7 /*SDL_GL_STENCIL_SIZE*/] > 0)
      };
      
	var ctx;// = Browser.createContext(canvas, is_SDL_OPENGL, true, webGLContextAttributes);
        var contextAttributes = {
          antialias: false,
          alpha: false
        };

        if (webGLContextAttributes) {
          for (var attribute in webGLContextAttributes) {
            contextAttributes[attribute] = webGLContextAttributes[attribute];
          }
        }

        var contextHandle = GL.createContext(canvas, contextAttributes);
        if (contextHandle) {
          ctx = GL.getContext(contextHandle).GLctx;
        }
        Module.ctx = ctx;
        GL.makeContextCurrent(contextHandle);
        Module.useWebGL = true;
        Browser.moduleContextCreatedCallbacks.forEach(function(callback) { callback() });
        Browser.init();
            
      SDL.surfaces[surf] = {
        width: width,
        height: height,
        canvas: canvas,
        ctx: ctx,
        surf: surf,
        buffer: buffer,
        pixelFormat: pixelFormat,
        alpha: 255,
        flags: flags,
        locked: 0,
        usePageCanvas: true,
        source: 'screen',

        isFlagSet: function(flag) {
          return flags & flag;
        }
      };

	SDL.screen = surf;
	
	return contextHandle;
    }
});
