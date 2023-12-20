// This produces adc, sp, sc for tests of DepoFluxWriter.  Note, this
// produces PER APA DATA PRODUCTS AND NOT MONOLITHS.  This is a more
// proper pattern for DUNE FD modules at least up to and including SP.

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/experiment/dune-vd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/dune10kt-1x2x6/simparams.jsonnet';
local fcl_params = {
    G4RefTime: std.extVar('G4RefTime') * wc.us,
};
local params = params_maker(fcl_params) {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
};


local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune10kt-1x2x6/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local output = 'wct-sim-ideal-sig.npz';
local sed_label = 'IonAndScint';

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  deposet: g.pnode({
    type: 'wclsSimDepoSetSource',
    name: "",
    data: {
      model: "",
      scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
      art_tag: sed_label, //name of upstream art producer of depos "label:instance:processName"
      assn_art_tag: "",
      id_is_track: false,     // use id-is-index trick
      debug_file: "wcls-sim-drift-deposource.log",
    },
  }, nin=0, nout=1),

};

// These type-names must be added to "outputers"
local save_digits = function(anode)
    local sident = '%d' % anode.data.ident;
    g.pipeline([
        g.pnode({
            type: 'wclsFrameSaver',
            name: 'adc%02d' % anode.data.ident,
            data: {
                anode: wc.tn(anode),
                digitize: true,
                frame_tags: ['orig' + sident],
                // nticks: params.daq.nticks,
                // chanmaskmaps: ['bad'],
                pedestal_mean: 'native',
            },
        }, nin=1, nout=1, uses=[anode]),
        g.pnode({
            type: 'DumpFrames',
            name: 'adc' + sident,
        }, nin=1, nout=0),
    ]);
local tap_digits(anode) = g.fan.tap('FrameFanout', save_digits(anode),
                                    name='adc%d'%anode.data.ident);


local save_signals = function(anode)
    local sident = '%d' % anode.data.ident;
    g.pipeline([
        g.pnode({
            type: 'wclsFrameSaver',
            name: 'sig%02d' % anode.data.ident,
            data: {
                anode: wc.tn(anode),
                digitize: false,
                frame_tags: ['gauss'+sident, 'wiener'+sident,'dnnsp'+sident],
                // frame_scale: [0.005, 0.005, 0.005], // utterly bogus thing to do
                chanmaskmaps: [],
                nticks: params.daq.nticks,
            },
        }, nin=1, nout=1, uses=[anode]),
        g.pnode({
            type: 'DumpFrames',
            name: 'sig' + sident,
        }, nin=1, nout=0)
    ]);

// save "threshold" from normal decon for each channel noise
// used in imaging
// sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),


local drifter = sim.drifter;
local setdrifter = g.pnode({
    type: 'DepoSetDrifter',
    data: {
        drifter: "Drifter"
    }
}, nin=1, nout=1, uses=[drifter]);

local depofluxwriter = g.pnode({
    type: 'wclsDepoFluxWriter',
    name: 'postdrift',          // must go into "outputers"
    data: {
        simchan_label: 'simpleSC',  // where to save in art::Event
        anodes: [wc.tn(anode) for anode in tools.anodes],
        field_response: wc.tn(tools.field),
        tick: 0.5 * wc.us,
        window_start: 0,
        window_duration: self.tick * 6000,
        nsigma: 3.0,
        // production usage should NOT set sed_label nor debug_file!
        sed_label: sed_label,
        debug_file: "wcls-sim-drift-depowriter.log",
    }
}, nin=1, nout=1, uses=tools.anodes+[tools.field]);


local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/dune10kt-1x2x6/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: '%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n),
  uses: [tools.anodes[n], tools.field],  // pnode extension
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/dune10kt-1x2x6/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/dune10kt-1x2x6/sp.jsonnet';
local sp = sp_maker(params, tools, { sparse: true });
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];



local multipass = [
    g.pipeline([
        // wcls_simchannel_sink[n],
        sn_pipes[n],
        tap_digits(tools.anodes[n]),
        // nf_pipes[n],
        sp_pipes[n],
        save_signals(tools.anodes[n]),
    ], 'multipass%d' % n)
    for n in anode_iota
];
local outtags = [];
local tag_rules = {
    frame: {
        '.*': 'framefanin',
    },
    trace: {['gauss%d' % anode.data.ident]: ['gauss%d' % anode.data.ident] for anode in tools.anodes}
        + {['wiener%d' % anode.data.ident]: ['wiener%d' % anode.data.ident] for anode in tools.anodes}
        + {['threshold%d' % anode.data.ident]: ['threshold%d' % anode.data.ident] for anode in tools.anodes}
        + {['dnnsp%d' % anode.data.ident]: ['dnnsp%d' % anode.data.ident] for anode in tools.anodes},
};




local fan = f.multifanout('DepoSetFanout', multipass,
                          fout_nnodes=[1,2], fout_multi=[2,6]);

// bi_manifold = f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,2], [2,6], [1,2], [2,6], 'sn_mag_nf', outtags, tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'gauss\\d+': 'gauss',
        'wiener\\d+': 'wiener',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
//local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.deposet, setdrifter, depofluxwriter, fan]);
//, bi_manifold, retagger, wcls_output.sp_signals, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
