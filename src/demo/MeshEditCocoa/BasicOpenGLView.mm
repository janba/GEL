#import "BasicOpenGLView.h"
#include <GEL/GLGraphics/MeshEditor.h>

using namespace CGLA;
using namespace GLGraphics;

MeshEditor me;

// ==================================
#pragma mark ---- Error Reporting ----

// error reporting as both window message and debugger string
void reportError (char * strError, int where)
{
    NSMutableDictionary *attribs = [NSMutableDictionary dictionary];
    [attribs setObject: [NSFont fontWithName: @"Monaco" size: 9.0f] forKey: NSFontAttributeName];
    [attribs setObject: [NSColor whiteColor] forKey: NSForegroundColorAttributeName];
    
    //	gErrorTime = getElapsedTime ();
	NSString * errString = [NSString stringWithFormat:@"Error: %s @ line %d", strError, where ];
	NSLog (@"%@\n", errString);
}

// ---------------------------------

// if error dump gl errors to debugger string, return error
GLenum glReportError (int where = -1)
{
	GLenum err = glGetError();
	if (GL_NO_ERROR != err)
		reportError ((char *) gluErrorString (err), where);
	return err;
}

#pragma mark ---- OpenGL Utils ----

// ---------------------------------


// ===================================

@implementation BasicOpenGLView

// ---------------------------------

// handles resizing of GL need context update and if the window dimensions change, a
// a window dimension update, reseting of viewport and an update of the projection matrix
- (void) reshape
{
    NSRect backRect = [self convertRectToBacking:[self bounds]];
    glViewport (0, 0, NSWidth(backRect), NSHeight(backRect));
	me.reshape(NSWidth(backRect),NSHeight(backRect));
}

// ---------------------------------
// ---------------------------------

// per-window timer function, basic time based animation preformed here
- (void)animationTimer:(NSTimer *)timer
{
    if(me.try_spinning_ball())
        [self setNeedsDisplay: YES];
}

-(IBAction)save_window_to_pasteboard:(id)sender
{
    // Get the pasteboard.
    NSPasteboard *pb = [NSPasteboard generalPasteboard];
    
    // Telling the pasteboard what type of data we're going to send in
    [pb declareTypes:[NSArray arrayWithObjects:NSPasteboardTypePNG,nil] owner:self];
    
    // Get the size of the image in a retina safe way
    NSRect backRect = [self convertRectToBacking: [self bounds]];
    int W = NSWidth(backRect);
    int H = NSHeight(backRect);

    // Create image. Note no alpha channel. I don't copy that.
    NSBitmapImageRep *rep = [[NSBitmapImageRep alloc] initWithBitmapDataPlanes: NULL
                                                                    pixelsWide: W
                                                                    pixelsHigh: H
                                                                 bitsPerSample: 8
                                                               samplesPerPixel: 3
                                                                      hasAlpha: NO
                                                                      isPlanar: NO
                                                                colorSpaceName: NSCalibratedRGBColorSpace
                                                                   bytesPerRow: 3*W
                                                                  bitsPerPixel: 0];
    
    // The following block does the actual reading of the image
    glPushAttrib(GL_PIXEL_MODE_BIT); // Save state about reading buffers
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1); // Dense packing
    glReadPixels(0, 0, W, H, GL_RGB, GL_UNSIGNED_BYTE, [rep bitmapData]);
    glPopAttrib();
    
    // So we need one more image, since we must flip its orientation.
    NSBitmapImageRep* flipped =  [[NSBitmapImageRep alloc] initWithBitmapDataPlanes: NULL
                                                                         pixelsWide: W
                                                                         pixelsHigh: H
                                                                      bitsPerSample: 8
                                                                    samplesPerPixel: 3
                                                                           hasAlpha: NO
                                                                           isPlanar: NO
                                                                     colorSpaceName: NSCalibratedRGBColorSpace
                                                                        bytesPerRow: 3*W
                                                                       bitsPerPixel: 0];
    
    // Primitive double for loop flipping the row order. Should be a better way. Can't find it.
    for(int j=1; j< H+1; ++j)
        for(int i=0;i<W;++i)
        {
            NSUInteger pixels[4];
            [rep getPixel: pixels atX:i y:j];
            [flipped setPixel: pixels atX:i y:H-j];
        }
    
    // Converting the representation to PNG and sending it to the pasteboard (with type indicated)
    [pb setData:[flipped representationUsingType:NSPNGFileType properties:nil] forType:NSPasteboardTypePNG];
}



-(void)keyDown:(NSEvent *)theEvent
{
    switch ([theEvent keyCode]) {
        case 123:
            me.key_left();
            break;
        case 124:
            me.key_right();
            break;
        case 126:
            me.key_up();
            break;
        case 125:
            me.key_down();
            break;
        case 115:
            me.key_home();
            break;
        case 119:
            me.key_end();
            break;
        default:
            NSString *characters = [theEvent characters];
            if ([characters length]) {
                unichar character = [characters characterAtIndex:0];
                me.keyparse(character);
            }
            break;
    }
    [self setNeedsDisplay: YES];
}


// ---------------------------------

- (void)mouseDown:(NSEvent *)theEvent // trackball
{
    if ([theEvent modifierFlags] & NSAlternateKeyMask) // send to pan
		[self rightMouseDown:theEvent];
	else {
        
        NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
        [[self openGLContext] makeCurrentContext];
        Vec2i pos(location.x,location.y);
        pos *= int([[self window]backingScaleFactor]);
        if([theEvent modifierFlags] & NSShiftKeyMask)
            me.grab_mesh(pos);
        else if([theEvent modifierFlags] & NSControlKeyMask) {
            if(me.select_vertex(pos))
                [self setNeedsDisplay: YES];
        }
        else
            me.grab_ball(ROTATE_ACTION,pos);
	}
}



// ---------------------------------

- (void)rightMouseDown:(NSEvent *)theEvent // pan
{
	NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
    Vec2i pos(location.x,location.y);
    pos *= int([[self window]backingScaleFactor]);
    me.grab_ball(PAN_ACTION,pos);
}

// ---------------------------------

- (void)otherMouseDown:(NSEvent *)theEvent //dolly
{
}

// ---------------------------------

- (void)mouseUp:(NSEvent *)theEvent
{
    me.release_mesh();
    me.release_ball();
}

// ---------------------------------

- (void)rightMouseUp:(NSEvent *)theEvent
{
	[self mouseUp:theEvent];
}

// ---------------------------------

- (void)otherMouseUp:(NSEvent *)theEvent
{
	[self mouseUp:theEvent];
}

// ---------------------------------

- (void)mouseDragged:(NSEvent *)theEvent
{
	NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
    Vec2i pos(location.x,location.y);
    pos *= int([[self window]backingScaleFactor]);
//    NSLog(@"Pos <%d %d> ", pos[0], pos[1]);
    [[self openGLContext] makeCurrentContext];
    if(!me.drag_mesh(pos))
        me.roll_ball(pos);
    [self setNeedsDisplay: YES];
    
}

// ---------------------------------

- (void)scrollWheel:(NSEvent *)theEvent
{
    Vec2i pos(0,300);
    me.grab_ball(ZOOM_ACTION,pos);
    pos[1] -= 3*[theEvent deltaY];
    me.roll_ball(pos);
    me.release_ball();
    [self setNeedsDisplay: YES];
}

// ---------------------------------

- (void)rightMouseDragged:(NSEvent *)theEvent
{
	[self mouseDragged: theEvent];
}

// ---------------------------------

- (void)otherMouseDragged:(NSEvent *)theEvent
{
	[self mouseDragged: theEvent];
}

// ---------------------------------

- (void) drawRect:(NSRect)rect
{
    NSOpenGLContext* ctxt = [self openGLContext];
    if(ctxt != nil) {
        [ctxt makeCurrentContext];
        me.display(int([[self window]backingScaleFactor]));
        [ctxt flushBuffer];
    }
}

// ---------------------------------

// set initial OpenGL state (current context is set)
// called after context is created
- (void) prepareOpenGL
{
    [[self openGLContext] makeCurrentContext];
    GLint swapInt = 1;
    [[self openGLContext] setValues:&swapInt forParameter:NSOpenGLCPSwapInterval]; // set to vbl sync
    me.init();
    NSLog(@"OpenGL Initialized");
}
// ---------------------------------

// this can be a troublesome call to do anything heavyweight, as it is called on window moves, resizes, and display config changes.  So be
// careful of doing too much here.
- (void) update // window resizes, moves and display changes (resize, depth and display config change)
{
	[super update];
	if (![self inLiveResize])  {// if not doing live resize
	}
    
}

// ---------------------------------

-(id) initWithFrame: (NSRect) frameRect
{
	self = [super initWithFrame: frameRect ];
    return self;
}

// ---------------------------------

- (BOOL)acceptsFirstResponder
{
    return YES;
}

// ---------------------------------

- (BOOL)becomeFirstResponder
{
    return  YES;
}

// ---------------------------------

- (BOOL)resignFirstResponder
{
    return YES;
}

// ---------------------------------

- (void) awakeFromNib
{
	timer = [NSTimer timerWithTimeInterval:(1.0f/60.0f) target:self selector:@selector(animationTimer:) userInfo:nil repeats:YES];
	[[NSRunLoop currentRunLoop] addTimer:timer forMode:NSDefaultRunLoopMode];
	[[NSRunLoop currentRunLoop] addTimer:timer forMode:NSEventTrackingRunLoopMode]; // ensure timer fires during resize
    [[NSApplication sharedApplication] setDelegate: self];
    [self registerForDraggedTypes: [NSArray arrayWithObjects: (id) kUTTypeData,nil]];
}

-(IBAction)open_file_dialog:(id)sender
{
    NSOpenPanel* openDlg = [NSOpenPanel openPanel];
    [openDlg setCanChooseFiles:YES];
    [openDlg setCanChooseDirectories:NO];
    [openDlg setAllowsMultipleSelection:YES];
    if ( [openDlg runModal] == NSOKButton )
        for(NSURL *fileURL in [openDlg URLs]) {
            me.add_file([ [fileURL path] UTF8String]);
            [[NSDocumentController sharedDocumentController] noteNewRecentDocumentURL: fileURL];
        }
    [self setNeedsDisplay: YES];
}

- (BOOL)application:(NSApplication *)theApplication openFile:(NSString *)filename
{
    if(me.add_file([ filename UTF8String])) {
        [[NSDocumentController sharedDocumentController] noteNewRecentDocumentURL:[NSURL fileURLWithPath:filename]];
        [self setNeedsDisplay: YES];
        return YES;
    }
    return NO;
}

- (BOOL)performDragOperation:(id <NSDraggingInfo>)sender
{
    NSURL* fileURL = [NSURL URLFromPasteboard: [sender draggingPasteboard]];
    me.add_file([[fileURL path] UTF8String]);
    [[NSDocumentController sharedDocumentController] noteNewRecentDocumentURL: fileURL];
    [self setNeedsDisplay: YES];
    return YES;
}

- (BOOL)prepareForDragOperation:(id <NSDraggingInfo>)sender
{
    return YES;
}

- (NSDragOperation)draggingEntered:(id <NSDraggingInfo>)sender
{
    return NSDragOperationCopy;
}


@end
