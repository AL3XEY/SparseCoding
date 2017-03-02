close all;
clear all;
if (is_octave)
  warning('off', 'Octave:missing-semicolon');
  warning('off', 'Octave:possible-matlab-short-circuit-operator');
  do_braindead_shortcircuit_evaluation (1);
  warning('off', 'Octave:associativity-change');
  warning('off', 'Octave:precedence-change');
  warning('off', 'Octave:str-to-num');
  warning('off', 'Octave:string-concat');
  
  warning("off", "all")
end
demo_fast_sc(1)