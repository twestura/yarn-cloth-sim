//
//  AppDelegate.mm
//  Simulator
//
//  Created by eschweickart on 2/10/14.
//
//

#import "AppDelegate.h"
#include "Simulator.h"

@implementation AppDelegate

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
  // Insert code here to initialize your application
  Simulator *s = new Simulator();
  s->run();
}

- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)AppDelegate
{
  return YES;
}

@end
