# EvolvingHP
The traits of CTRNN-based homeostatic plasticity mechanisms are evolved to best maintain the function of a circuit approximating the pyloric rhythm.

Workflow is generally that several "main____.cpp" files are compiled with their corresponding "Makefile_____"'s and data is output to a file, so it can be imported to visualize in Jupyter notebooks or Mathematica notebooks.

Flow of experiment is planned as follows:

<ol>
  <li>Evolve 3-N CTRNN pyloric solution(s)</li>
  <li>Map parameter space in the vicinity of best solutions, then identify and explain distinct fitness regimes and boundaries.</li>
  <li>Pick paradigmatic initial point in each of the different regimes</li>
  <li>Evolve on an HP mechanism which “recovers” function best from the greatest number of initial conditions (HP mechanism Fitness = avg final pyloric fitness between all initial conditions)</li> 
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
    <li>How well do they align with fitness landscapes mapped from fast dynamics alone? (i.e. is “realm of acceptability” enough to explain success?)</li>
    <ol>
      <li>There very well could be parts where the HP itself has produced a fit solution that is not otherwise present, but this isn’t the main story for this arc</li>
    </ol>
    <li>Finite state machine of subsequent perturbations - even if you can’t recover from one, a subsequent perturbation might lead to full recovery</li>
    <ol>
      <li>In conversation with Connor’s viability mapping procedures</li>
    </ol>
  </ol>
</ol>
