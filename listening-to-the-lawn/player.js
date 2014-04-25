//adapted from jasmid Replayer
function Player(midiFile, fns) {
	var noteOn = fns.noteOn;
	var noteOff = fns.noteOff;
	var programChange = fns.programChange;
	var lastEvent = fns.lastEvent;
	if(noteOn==null) noteOn = function(){};
	if(noteOff==null) noteOff = function(){};
	if(programChange==null) programChange = function(){};
	if(lastEvent==null) lastEvent = function(){};
	
	var trackStates = [];
	var beatsPerMinute = 120;
	var ticksPerBeat = midiFile.header.ticksPerBeat;
	var channelCount = 16;
	var time = 0;
	
	for (var i = 0; i < midiFile.tracks.length; i++) {
		trackStates[i] = {
			'nextEventIndex': 0,
			'ticksToNextEvent': (
				midiFile.tracks[i].length ?
					midiFile.tracks[i][0].deltaTime :
					null
			)
		};
	}
		
	var nextEventInfo;
	var secondsToNextEvent = 0;
	
	function getNextEvent() {
		var ticksToNextEvent = null;
		var nextEventTrack = null;
		var nextEventIndex = null;
		
		for (var i = 0; i < trackStates.length; i++) {
			if (
				trackStates[i].ticksToNextEvent != null
				&& (ticksToNextEvent == null || trackStates[i].ticksToNextEvent < ticksToNextEvent)
			) {
				ticksToNextEvent = trackStates[i].ticksToNextEvent;
				nextEventTrack = i;
				nextEventIndex = trackStates[i].nextEventIndex;
			}
		}
		if (nextEventTrack != null) {
			/* consume event from that track */
			var nextEvent = midiFile.tracks[nextEventTrack][nextEventIndex];
			if (midiFile.tracks[nextEventTrack][nextEventIndex + 1]) {
				trackStates[nextEventTrack].ticksToNextEvent += midiFile.tracks[nextEventTrack][nextEventIndex + 1].deltaTime;
			} else {
				trackStates[nextEventTrack].ticksToNextEvent = null;
			}
			trackStates[nextEventTrack].nextEventIndex += 1;
			/* advance timings on all tracks by ticksToNextEvent */
			for (var i = 0; i < trackStates.length; i++) {
				if (trackStates[i].ticksToNextEvent != null) {
					trackStates[i].ticksToNextEvent -= ticksToNextEvent
				}
			}
			nextEventInfo = {
				'ticksToEvent': ticksToNextEvent,
				'event': nextEvent,
				'track': nextEventTrack
			}
			var beatsToNextEvent = ticksToNextEvent / ticksPerBeat;
			secondsToNextEvent = beatsToNextEvent / (beatsPerMinute / 60);
			time+=secondsToNextEvent;
		} else {
			nextEventInfo = null;
			self.finished = true;
		}
	}
	
	function next(){
		getNextEvent();	
		if(nextEventInfo!=null) {
			setTimeout(handleEvent, secondsToNextEvent*1000);
		}
		else{
			lastEvent();
		}
	}
	
	function play(samples) {	
		next();
	}
	
	function handleEvent() {
		var event = nextEventInfo.event;
		switch (event.type) {
			case 'meta':
				switch (event.subtype) {
					case 'setTempo':
						beatsPerMinute = 60000000 / event.microsecondsPerBeat
				}
				break;
			case 'channel':
				switch (event.subtype) {
					case 'noteOn':
						noteOn(time, event.channel, event.noteNumber, event.velocity);
						break;
					case 'noteOff':
						noteOff(time, event.channel, event.noteNumber, event.velocity);
						break;
					case 'programChange':
						programChange(time, event.programNumber);
				}
				break;
		}
		next();
	}

	var self = {
		'play': play,
		'finished': false
	}
	return self;
}