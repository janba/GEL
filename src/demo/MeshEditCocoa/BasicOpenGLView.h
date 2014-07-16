@interface BasicOpenGLView : NSOpenGLView < NSDraggingDestination, NSApplicationDelegate >
{
	NSTimer* timer;
 
    bool fAnimate;
	IBOutlet NSMenuItem * animateMenuItem;
    bool fInfo;
	IBOutlet NSMenuItem * infoMenuItem;
	
	CFAbsoluteTime time;
}

-(IBAction)save_window_to_pasteboard:(id)sender;

- (void) keyDown:(NSEvent *)theEvent;
- (void) mouseDown:(NSEvent *)theEvent;
- (void) rightMouseDown:(NSEvent *)theEvent;
- (void) otherMouseDown:(NSEvent *)theEvent;
- (void) mouseUp:(NSEvent *)theEvent;
- (void) rightMouseUp:(NSEvent *)theEvent;
- (void) otherMouseUp:(NSEvent *)theEvent;
- (void) mouseDragged:(NSEvent *)theEvent;
- (void) scrollWheel:(NSEvent *)theEvent;
- (void) rightMouseDragged:(NSEvent *)theEvent;
- (void) otherMouseDragged:(NSEvent *)theEvent;

- (void) drawRect:(NSRect)rect;

-(IBAction)open_file_dialog:(id)sender;

- (void) prepareOpenGL;
- (void) update;		// moved or resized

- (BOOL) acceptsFirstResponder;
- (BOOL) becomeFirstResponder;
- (BOOL) resignFirstResponder;

- (id) initWithFrame: (NSRect) frameRect;
- (void) awakeFromNib;
- (BOOL)performDragOperation:(id <NSDraggingInfo>)sender;
@end
