//
//  main.m
//  ImageTest
//
//  Created by Simon Frost on 30/06/2014.
//  Copyright (c) 2014 Orangeninja. All rights reserved.
//

#import <Cocoa/Cocoa.h>

int main(int argc, const char * argv[])
{
  return NSApplicationMain(argc, argv);
  
  double pixelRatio = 2.4;
  
  for (int col=0; col<4; col++) {
    double realX = (double)col * pixelRatio;
    int startX = floor(realX);
    int endX = floor(realX + pixelRatio);
    
    NSLog(@"Processing pixel %d, composed of pixels from column %d to %d", col, startX, endX);
    
    double total = 0;
    for (int srcX=startX; srcX<=endX; srcX++) {
      double frac;
      
      frac = 1;
      if (srcX == startX) {
        frac = fmod(srcX+1, pixelRatio);
      } else if (srcX == endX) {
        frac = 1-fmod(srcX+1, pixelRatio);
      }
      
      
      
      
//      double frac = srcX+1 % pixelRatio;
/*      double frac = fmod(srcX, pixelRatio);
      if (srcX > realX + pixelRatio) {
        frac = 1 - frac;
      }
      frac = MIN(frac, 1.0);
*/

      
      /*
      frac = fmod(srcX+1, pixelRatio);
      if (srcX+1 > realX + pixelRatio) {
        frac = 1 - frac;
      }
      frac = MIN(frac, 1.0);
      */
      
      
      
      NSLog(@"  Pixel %d contribution: %f", srcX, frac);
      total += frac;
    }
    NSLog(@"  Total: %f", total);
    
    NSLog(@" ");
  }
}
