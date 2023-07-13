st=>start: Start
e=>end
io1=>inputoutput: input: cls, stream, window_s, windowtype, inv, data_cleaner, z_threshold
op1=>operation: Align stream traces
op2=>operation: Adjust window samples to next highest power of 2
op2a=>operation: Calculate spectral multiplication factor
mulfact = 2 / (window_samples * sample_rate)
op3=>operation: Set up arrays to hold Fourier transforms, instrument responses, and units
op4=>operation: Make a list of trace ids and select the first one
op6=>operation: Increment the trace id
op7=>operation: Select the corresponding trace
sub1=>subroutine: _calculate_windowed_rfft()
sub2=>subroutine: _correct_response()
cond1=>condition: At last trace?
sub3=>subroutine: _remove_outliers()
cond2=>condition: was a DataCleaner provided?
op8=>operation: Apply data cleaner, update channel names, update response channel keys
op9=>operation: Create a SpectralDensity object
op10=>operation: Set input and output channels to first channel
op11=>operation: Calculate cross-spectral density
x = 2*mean(conj(ft(input))*ft(output))*mulfact
mean is across windows
op12=>operation: put xsd in its SpectralDensity bin
cond3=>condition: At last output channel?
cond4=>condition: At last input channel?
op13=>operation: increment output channel
op14=>operation: increment input channel

st->io1->op1->op2->op2a->op3->op4
op4->op7->sub1->sub2->cond1
cond1(yes)->sub3->cond2
cond1(no)->op6->op7
cond2(yes, right)->op8->op9
cond2(no)->op9
op9->op10->op11->cond3
cond3(no)->op13->op11
cond3(yes)->cond4
cond4(no)->op14->op11
cond4(yes)->e
