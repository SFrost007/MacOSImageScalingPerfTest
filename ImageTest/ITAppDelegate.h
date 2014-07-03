//
//  ITAppDelegate.h
//  ImageTest
//
//  Created by Simon Frost on 30/06/2014.
//  Copyright (c) 2014 Orangeninja. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface ITAppDelegate : NSObject <NSApplicationDelegate, NSWindowDelegate>

@property (assign) IBOutlet NSWindow *window;
@property (weak) IBOutlet NSSegmentedControl *algorithmSwitch;
@property (weak) IBOutlet NSImageView *imageView;
@property (weak) IBOutlet NSTextField *computeTimeLabel;
@property (weak) IBOutlet NSTextField *resolutionLabel;

- (IBAction)algorithmChanged:(id)sender;

@end
