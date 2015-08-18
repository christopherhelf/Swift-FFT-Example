//
//  main.swift
//  ffttest
//
//  Created by Christopher Helf on 17.08.15.
//  Copyright (c) 2015 Christopher Helf. All rights reserved.
//

import Foundation
import Accelerate

var fft = FFT()

let n = 512 // Should be power of two for the FFT
let frequency1 = 4.0
let phase1 = 0.0
let amplitude1 = 8.0
let seconds = 2.0
let fps = Double(n)/seconds

var sineWave = (0..<n).map {
    amplitude1 * sin(2.0 * M_PI / fps * Double($0) * frequency1 + phase1)
}

fft.calculate(sineWave, fps: fps)

