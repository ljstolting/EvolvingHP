# EvolvingHP
The traits of CTRNN-based homeostatic plasticity mechanisms are evolved to best maintain the function of a circuit approximating the pyloric rhythm.

Workflow is generally that several "main____.cpp" files are compiled with their corresponding "Makefile_____"'s and data is output to a file, so it can be imported to visualize in Jupyter notebooks or Mathematica notebooks.

An HP genome will begin by specifying which of the 12 possible weight and bias parameters (in a 3N CTRNN) are subject to change via the HP mechanism (Williams, 2005). Then, it specifies the time-constants, lower bounds, ranges, and sliding windows, of each relevant neuron (in that order). In the future, it might become necessary to specify these separately for each parameter, but this is parsimonious enough for now and seems to be effective. It is equivalent to the assumption that there is only one "detector" mechanism for abberant activity in a neuron, and that this signal is simply equivalently fed to mechanisms of change for intramembrane ion channels and synaptic proteins (though synaptic weights still change proporitonately to their current value, as in Williams 2005 -- should I change that?). The structure of an HP genome is thus as follows:

bool_bias1 bool_bias2 bool_bias3 bool_w11 bool_w12 ... bool_w33

{genotype: all three neurons have at least one attribute changing}
HPtau_n1 HPtau_n2 HPtau_n3 lower_n1 lower_n2 lower_n3 range_n1 range_n2 range_n3 swindow_n1 swindow_n2 swindow_n3
{alt genotype: only neurons 1 & 3 have at least one attribute changing}
HPtau_n1 HPtau_n3 lower_n1 lower_n3 range_n1 range_n3 swindow_n1 swindow_n3

{phenotype}
""

Note that the neuron identities are as follows: 1-LP, 2-PY, 3-PD

Flow of inquiry is planned as follows:

<ol>
  <li>Evolve 3-N CTRNN pyloric solution(s)</li>
  <li>Map parameter space in the vicinity of best solutions, then identify and explain distinct fitness regimes and boundaries.</li>
  <li>Pick paradigmatic initial points in each of the different regimes around a successful solution</li>
  <li>Evolve on an HP mechanism which “recovers” function best from the greatest number of initial conditions (HP mechanism Fitness = avg final pyloric fitness between all initial conditions) Note that pyloricness must be evaluated at several time points to be sure that the solutions stay in the pyloric region, and are not dependent on the specific amount of time HP is given to activate</li> 

  Future Questions:
  <ol>
    <li>Are there certain ones that simply can’t be recovered/are harder to recover?</li>
    <li>Are there certain tradeoffs between different types of recovery proficiency</li>
    <ol>
      <li>Might organisms choose to be more reliable in cases they expect to experience more often?</li>
    </ol>
    <li>How can you explain the degree of success based on what information is available to the HP mechanism</li>
    <ol>
      <li>“Capital I” Information about behavioral fitness based on average activity level (which is not perfectly informative)</li>
    </ol>
    <li>Does it generalize to other similar types of activity profiles, in different (even discontinuous) regions of parameter space?</li>
    <li>Do the best HP mechanisms generalize to other (degenerate) pyloric solutions?</li>
    <ol>
      <li>Would generalization be better when the way I define behavioral fitness isn’t so arbitrary and discrete? (i.e. with the walker)</li>
    </ol>
  </ol>
  <li>For a given evolved mechanism, map stabilities in higher-dimensional parameter space (probably mostly limit cycles, but occasionally might approximate static)</li>
  <ol>
    <li>Are the timescales of HP separate enough to explain success by the "realm of acceptability" paradigm (where on average the push to leave that point is low)</li>
    <ol>
      <li>There very well could be parts where the HP itself has produced a fit solution that is not otherwise present, but this isn’t the main story for this arc</li>
    </ol>
    <li>Finite state machine of subsequent perturbations - even if you can’t recover from one, a subsequent perturbation might lead to full recovery</li>
    <ol>
      <li>In conversation with Connor’s viability mapping procedures</li>
    </ol>
  </ol>
</ol>
