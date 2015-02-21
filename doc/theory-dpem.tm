<TeXmacs|1.99.2>

<style|generic>

<\body>
  <doc-data|<doc-title|Theoretical Justification for the EM-Style Inference
  of DP Mixture Models>|<doc-author|<author-data|<author-name|Richard
  Guo>|<\author-affiliation>
    Department of Computer Science, Duke University
  </author-affiliation>|<\author-affiliation>
    <samp|guo@cs.duke.edu>
  </author-affiliation>>>>

  <section|Collapsing CRP into finite-dimensional <math|\<pi\><rsub|K>>>

  In this section, we introduce a finite dimensional partition distribution
  <math|\<pi\><rsub|K>>, which is constructed by keeping the largest <math|K>
  clusters in the partition distribution induced by DP (i.e. the CRP).\ 

  There are two representations of DP mixture models, one through the
  distribution of sample-specific parameters, one through the CRP partition
  distribution. The former representation is usually written as\ 

  <\equation*>
    \<Theta\>\<sim\>DP<around*|(|\<alpha\>,G<rsub|0>|)><space|1em>\<theta\><rsub|1>,\<ldots\>,\<theta\><rsub|n><above|\<sim\>|iid>\<Theta\><space|1em>X<rsub|i>\|\<theta\><rsub|i><above|\<sim\>|ind>P<around*|(|X\|\<theta\><rsub|i>|)>.
  </equation*>

  Because <math|\<Theta\>> is an atomic distribution almost surely, which
  allows the representation <math|\<Theta\>=<big|sum><rsub|k=1><rsup|\<infty\>>\<gamma\><rsub|k>
  \<delta\><rsub|\<theta\><rsup|\<ast\>><rsub|k>>>, the samples
  <math|<around*|{|\<theta\><rsub|i>|}>> usually have repetitions, i.e.
  clusters. The second representation is a more straightforward description
  of this perspective. Let us say each <math|X<rsub|i>> is associated with a
  latent variable <math|Z<rsub|i>\<in\>\<bbb-Z\><rsup|+>> that numbers its
  cluster. They are drawn from the partition distribution induced from DP,
  which is also referred to as the CRP (Chinese Restaurant Process)\ 

  <\equation*>
    <around*|{|Z<rsub|1>,\<ldots\>,Z<rsub|n>|}>\<sim\>CRP<around*|(|\<alpha\>|)>\<nocomma\>.
  </equation*>

  Note that these draws are not independent because otherwise CRP cannot
  induce the desired ``rich-get-richer'' clustering effect. Fixing <math|n>,
  The partition can be expressed as a function
  <math|f:<around*|{|1,\<ldots\>,n|}>\<rightarrow\>\<bbb-Z\><rsub|+>> that
  maps each <math|Z<rsub|i>> to a natural number. For a finite <math|n>,
  <math|CRP<around*|(|\<cdummy\>|)>> can be viewed as a distribution over all
  such <math|f>. This representation decouples the clustering structure and
  the drawing of cluster-specific parameters by

  <\equation*>
    <around*|{|\<theta\><rsub|c>|}> <above|\<sim\>|iid>
    G<rsub|0>,<space|1em>P<around*|(|X<rsub|i>\|<around*|{|Z<rsub|j>|}>|)>=P<around*|(|X<rsub|i>\|Z<rsub|i>|)>=P<around*|(|X<rsub|i>\|\<theta\><rsub|Z<rsub|i>>|)>.
  </equation*>

  For such a draw <math|<around*|{|Z<rsub|i>|}><rsub|i=1><rsup|n>>, imagine
  we explicitly sort the size of clusters in decreasing order and number the
  biggest cluster as 1, the 2nd biggest as 2, etc. Then
  <math|<around*|{|Z<rsub|1>,\<ldots\>,Z<rsub|n>|}>> can take the values
  contiguously from 1 to at most <math|n>. Given an integer
  <math|K\<geqslant\>1>, let us keep the first <math|K> clusters in
  <math|<around*|{|Z<rsub|i>|}>> and meanwhile merging all the other smaller
  clusters into <with|font-series|bold|a single cluster> called
  <math|\<ast\>>, as shown in Figure <reference|fig:pik>. The resulting
  partition is a mapping <math|g<rsub|K>:<around*|{|1,\<ldots\>,n|}>\<rightarrow\><around*|{|1,\<ldots\>,K|}>\<cup\><around*|{|\<ast\>|}>>.
  We call the <with|font-series|bold|induced> distribution on
  <math|g<rsub|K>> as <math|\<pi\><rsub|K>>, a
  <math|<around*|(|K+1|)>>-dimensional partition distribution over <math|n>
  samples. Notice that (1) empty clusters are allowed in a
  <math|\<pi\><rsub|K>> draw (2) the order among clusters in a CRP draw can
  be set arbitrarily.

  <big-figure|<with|gr-mode|<tuple|edit|text-at>|gr-frame|<tuple|scale|1cm|<tuple|0.5gw|0.519999gh>>|gr-geometry|<tuple|geometry|1par|0.6par>|gr-grid|<tuple|empty>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|empty>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-auto-crop|true|gr-crop-padding|0spc|<graphics||<point|-6|0>|<point|-5.5|0>|<point|-5|0>|<point|-4.4|0>|<point|-4|0>|<point|-3.5|0>|<point|-3|0>|<point|-1.3|0>|<point|-0.8|0>|<point|-0.3|0>|<point|0.2|0>|<point|3|0>|<point|3.4|0>|<point|4.7|0>|<point|5.7|0>|<point|7|0>|<line|<point|-6.4|0.5>|<point|-2.5|0.5>|<point|-2.5|-0.5>|<point|-6.4|-0.5>|<point|-6.4|0.5>>|<line|<point|-1.6|0.4>|<point|0.599999999999999|0.4>|<point|0.599999999999999|-0.5>|<point|-1.6|-0.5>|<point|-1.6|0.4>>|<line|<point|2.7|0.3>|<point|3.7|0.3>|<point|3.7|-0.3>|<point|2.7|-0.3>|<point|2.7|0.3>>|<line|<point|4.4|0.3>|<point|5.0|0.3>|<point|5.0|-0.3>|<point|4.4|-0.3>|<point|4.4|0.3>>|<line|<point|5.4|0.3>|<point|6.0|0.3>|<point|6.0|-0.3>|<point|5.4|-0.3>|<point|5.4|0.3>>|<point|-6|-1.4>|<point|-5.5|-1.4>|<point|-5|-1.4>|<point|-4.5|-1.4>|<point|-4|-1.4>|<point|-3.5|-1.4>|<point|-3|-1.4>|<point|-1.3|-1.4>|<point|-0.8|-1.4>|<point|-0.3|-1.4>|<point|0.3|-1.4>|<point|3|-1.4>|<point|3.4|-1.4>|<point|4.7|-1.4>|<point|5.7|-1.4>|<point|7|-1.4>|<line|<point|-6.3|-1>|<point|-2.5|-1.0>|<point|-2.5|-1.7>|<point|-6.3|-1.7>|<point|-6.3|-1.0>>|<line|<point|-1.5|-1>|<point|0.599999999999999|-1.0>|<point|0.599999999999999|-1.7>|<point|-1.5|-1.7>|<point|-1.5|-1.0>>|<line|<point|2.7|-1>|<point|7.3|-1.0>|<point|7.3|-1.8>|<point|2.7|-1.8>|<point|2.7|-1.0>>|<line|<point|6.7|0.3>|<point|7.3|0.3>|<point|7.3|-0.3>|<point|6.7|-0.3>|<point|6.7|0.3>>|<text-at|1|<point|-4.5|0.8>>|<text-at|2|<point|-0.6|0.8>>|<text-at|3|<point|3|0.7>>|<text-at|4|<point|4.6|0.7>>|<text-at|5|<point|5.6|0.7>>|<text-at|6|<point|7.0|0.7>>|<text-at|1|<point|-4.6|-2.2>>|<text-at|2|<point|-0.5|-2.2>>|<text-at|\<ast\>|<point|5|-2.2>>|<text-at|CRP|<point|-7.5|0>>|<text-at|<math|\<pi\><rsub|2>>|<point|-7.3|-1.6>>|<text-at||<point|1.51308|1.26511>>>>|An
  example where a draw from CRP is collapsed into a <math|\<pi\><rsub|2>>
  draw. <label|fig:pik>>

  Intuitively, <math|\<pi\><rsub|K>> is a filtered version of CRP which we
  only focus on the largest <math|K> clusters, throwing all the rest into a
  reservoir denoted by <math|\<ast\>>. Meanwhile, CRP can be viewed as the
  limiting distribution of <math|\<pi\><rsub|K>> when
  <math|K\<rightarrow\>\<infty\>>. In fact, given finite sample size
  <math|n>, <math|\<pi\><rsub|n>> is equivalent to CRP.\ 

  Directly characterizing CRP is hard, the likelihood for a particular
  partition involves an exponential number of additions and is largely
  intractable. However, arguably the most useful characterization is its full
  conditional distribution, which is also called the ``prediction rule'' by
  Ishwaran & James, namely

  <\equation*>
    P<around*|(|Z<rsub|i>\|<around*|{|Z<rsub|-i>|}>|)>=<frac|<big|sum><rsub|j\<neq\>i>\<delta\><rsub|Z<rsub|j>>+\<alpha\>\<delta\><rsub|\<ast\>>|n-1+\<alpha\>>=<frac|<big|sum><rsub|k>n<rsub|-i,k>
    \<delta\><rsub|k>+\<alpha\>\<delta\><rsub|\<ast\>>|n-1+\<alpha\>>,
  </equation*>

  where <math|\<ast\>> refers a new cluster that is different from all the
  existing clusters taken by <math|<around*|{|Z<rsub|-i>|}>>, and
  <math|<around*|{|Z<rsub|-i>|}>\<assign\><around*|{|Z<rsub|j>:1\<leqslant\>j\<leqslant\>n\<nocomma\>,j\<neq\>i|}>>.\ 

  Now we show that a similar property is inherited by <math|\<pi\><rsub|K>>,
  which is used for characterizing <math|\<pi\><rsub|K>> in our inference
  algorithm. Firstly, clusters <math|1\<leqslant\>k\<leqslant\>K> are
  one-to-one mapped between <math|\<pi\><rsub|K>> and CRP, and thus\ 

  <\equation*>
    \<pi\><rsub|K><around*|(|Z<rsub|i>=k\|<around*|{|Z<rsub|-i>|}>|)>=<frac|<big|sum><rsub|j\<neq\>i><with|font-series|bold|1><around*|(|Z<rsub|j>=k|)>|n-1+\<alpha\>><space|1em><around*|(|1\<leqslant\>k\<leqslant\>K|)>.
  </equation*>

  Secondly,\ 

  <\equation*>
    \<pi\><rsub|K><around*|(|Z<rsub|i>=\<ast\>\|<around*|{|Z<rsub|-i>|}>|)>=CRP<around*|(|Z<rsub|i>\<nin\><around*|{|1,\<ldots\>,K|}>\|<around*|{|Z<rsub|-i>|}>|)>=<frac|<around*|(|n-1|)>-<big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>\<in\><around*|{|1,\<ldots\>,K|}>|)>|n-1+\<alpha\>>+<frac|\<alpha\>|n-1+\<alpha\>>=<frac|n<rsub|-i,\<ast\>>+\<alpha\>|n-1+\<alpha\>>.
  </equation*>

  Therefore, the ``prediction rule'' for <math|\<pi\><rsub|K>> can be written
  as\ 

  \;

  <\equation*>
    \<pi\><rsub|K><around*|(|Z<rsub|i>\|<around*|{|Z<rsub|-i>|}>|)>=<frac|<big|sum><rsub|1\<leqslant\>k\<leqslant\>K>n<rsub|-i,k>
    \<delta\><rsub|k>+<around*|(|n<rsub|-i,\<ast\>>+\<alpha\>|)>\<delta\><rsub|\<ast\>>|n-1+\<alpha\>>.<label|eqs:pik-prediction>
  </equation*>

  Clearly, when <math|K\<geqslant\>n>, we have <math|n<rsub|-i,\<ast\>>=0>
  and the above becomes identical to CRP.

  <section|Overview of Objectives>

  The algorithm is structured as\ 

  <\render-code>
    Outer Loop: Increments <math|K> until <math|\<ast\>> is small enough

    {

    \ \ \ \ Inner Loop: Estimate the mixture model under
    <math|\<pi\><rsub|K>> with EM

    }
  </render-code>

  The objective is to show:

  <\enumerate-numeric>
    <item>The algorithm terminates with finite number of steps,

    <item>Each inner loop improves the marginal likelihood
    <math|P<around*|(|<around*|{|X<rsub|i>|}>\|\<pi\><rsub|K>|)>>,

    <item>In the outer loop, replacing <math|P<around*|(|<around*|{|X<rsub|i>|}>\|\<pi\><rsub|K>|)>>
    with <math|P<around*|(|<around*|{|X<rsub|i>|}>\|\<pi\><rsub|K+1>|)>> also
    improves that maximum marginal likelihood.<marginal-note|normal|c|For 3,
    can we argue by observing that the <math|\<pi\><rsub|K+1>>-mixture model
    space subsumes the <math|\<pi\><rsub|K>>-mixture?>
  </enumerate-numeric>

  <section|Theory for the EM procedure>

  The EM inner loop solves the parameter estimation problem for a
  <math|\<pi\><rsub|K>>-mixture model. Supposing we have data
  <math|<around*|{|X<rsub|i>|}><rsub|i=1><rsup|n>> and latent variables
  <math|<around*|{|Z<rsub|i>|}><rsub|i=1><rsup|n>> corresponding to each
  observation. Given the cluster label, the likelihood is available as\ 

  <\equation*>
    L<rsub|c><around*|(|X<rsub|i>|)>\<assign\>P<around*|(|X<rsub|i>\|Z<rsub|i>=c|)>=<choice|<tformat|<table|<row|<cell|P<around*|(|X<rsub|i>\|\<theta\><rsub|c>|)><space|1em><around*|(|1\<leqslant\>c\<leqslant\>K|)>>>|<row|<cell|<big|int>P<around*|(|X\|\<theta\><rsup|\<ast\>>|)>G<rsub|0><around*|(|\<theta\><rsup|\<ast\>>|)>d\<theta\><rsup|\<ast\>><space|1em><around*|(|c=\<ast\>|)>>>>>>.
  </equation*>

  The objective is to estimate\ 

  <\enumerate-roman>
    <item>the parameters <math|<around*|{|\<theta\><rsub|c>|}><rsub|c=1><rsup|K>>;

    <item>the distribution of latent cluster indicators
    <math|<around*|{|Z<rsub|i>|}>>. Yet, as the joint distribution is
    intractable and to some extent unnecessary, we only attempt the estimate
    the marginals <math|P<around*|(|Z<rsub|i>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>>.
  </enumerate-roman>

  The finite-dimensional mixture is the joint distribution
  <math|\<pi\><rsub|K><around*|(|<around*|{|Z<rsub|i>|}><rsub|i=1><rsup|n>|)>>,
  whereas each <math|z<rsub|i>\<in\><around*|{|1,\<ldots\>,K|}>\<cup\><around*|{|\<ast\>|}>>
  for a fixed number of components <math|K>.
  <math|\<pi\><around*|(|\<cdot\>|)>> is a finite-dimensional surrogate for
  the partition distribution induced by DP. Although a constructive
  characterization of <math|\<pi\><rsub|K><around*|(|\<cdot\>|)>> seems
  intractable, we know <math|\<pi\><rsub|K><around*|(|\<cdot\>|)>> should
  satisfy ``the prediction rule'' of CRP, i.e. the full conditional property
  of\ 

  <\equation*>
    \<pi\><rsub|K><around*|(|Z<rsub|i>\|<around*|{|Z<rsub|-i>|}>|)>=<frac|<big|sum><rsub|j\<neq\>i>\<delta\><rsub|Z<rsub|j>><around*|(|\<cdot\>|)>+<around*|(|\<alpha\>+n<rsub|-i,\<ast\>>|)>\<delta\><rsub|\<ast\>><around*|(|\<cdot\>|)>|n-1+\<alpha\>>=<choice|<tformat|<table|<row|<cell|<frac|<big|sum><rsub|j\<neq\>i><with|font-series|bold|1><around*|(|Z<rsub|j>=c|)>|n-1+\<alpha\>><space|1em><around*|(|1\<leqslant\>Z<rsub|i>\<leqslant\>K|)>>>|<row|<cell|<frac|\<alpha\>+<big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=\<ast\>|)>|n-1+\<alpha\>><space|1em><around*|(|Z<rsub|i>=\<ast\>|)>>>>>>.
  </equation*>

  By the usual logic of EM (or variational inference), introducing a
  distribution <math|q<around*|(|\<cdot\>|)>> over latent variables, the
  likelihood can be decomposed as

  <\eqnarray*>
    <tformat|<table|<row|<cell|log P<around*|(|<around*|{|X<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>>|<cell|=>|<cell|<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
    log<frac|P<around*|(|<around*|{|X<rsub|i>|}>,<around*|{|Z<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>|q<around*|(|<around*|{|Z<rsub|i>|}>|)>>+<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
    log<frac|q<around*|(|<around*|{|Z<rsub|i>|}>|)>|P<around*|(|<around*|{|Z<rsub|i>|}>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>><eq-number>>>|<row|<cell|>|<cell|=>|<cell|\<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>+H<around*|(|q|)>+KL<around*|(|q<around*|\|||\|>P<around*|(|<around*|{|Z<rsub|i>|}>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>|)>,<label|eqs:EM>>>>>
  </eqnarray*>

  where\ 

  <\equation*>
    \<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>=<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
    log P<around*|(|<around*|{|X<rsub|i>|}>,<around*|{|Z<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>
  </equation*>

  is the expected full-data log-likelihood under the
  <math|q<around*|(|\<cdot\>|)>>, <math|H<around*|(|q|)>=-<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
  log q<around*|(|<around*|{|Z<rsub|i>|}>|)>> is the entropy of
  <math|q<around*|(|\<cdot\>|)>> and the KL term is the difference between
  the introduced <math|q<around*|(|\<cdot\>|)>> and true posterior
  distribution of <math|<around*|{|Z<rsub|i>|}>>.\ 

  <subsection|The M-Step>

  We show that the M-step improves <math|\<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>>.
  <marginal-note|normal|c|Can also show the KL term is non-decreasing in the
  M-step?>

  Let us denote the marginal distribution of
  <math|q<around*|(|<around*|{|Z<rsub|i>|}>|)>> as
  <math|<around*|{|w<rsub|i,c>|}>>, i.e.\ 

  <\equation*>
    w<rsub|i,c>\<assign\><big|sum><rsub|<around*|{|Z<rsub|-i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>\|<rsub|Z<rsub|i>=c><space|1em><around*|(|c\<in\><around*|{|1,\<ldots\>,K|}>\<cup\><around*|{|\<ast\>|}>|)>
  </equation*>

  Remembering <math|P<around*|(|<around*|{|X<rsub|i>|}>,<around*|{|Z<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>=\<pi\><rsub|K><around*|(|<around*|{|Z<rsub|i>|}>|)>P<around*|(|<around*|{|X<rsub|i>|}>\|<around*|{|Z<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>>,
  we have \ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>>|<cell|=>|<cell|<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)><around*|(|log
    \<pi\><rsub|K><around*|(|<around*|{|Z<rsub|i>|}>|)>+log
    P<around*|(|<around*|{|X<rsub|i>|}>\|<around*|{|Z<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>log
    \<pi\><rsub|K><around*|(|<around*|{|Z<rsub|i>|}>|)>
    \<noplus\>\<noplus\>+<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)><big|sum><rsub|i=1><rsup|n><big|sum><rsub|c><with|font-series|bold|1><around*|(|Z<rsub|i>=c|)>
    log L<rsub|c><around*|(|X<rsub|i>\|\<theta\><rsub|c>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
    log \<pi\><rsub|K><around*|(|<around*|{|Z<rsub|i>|}>|)>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|c\<in\><around*|[|K|]>\<cup\><around*|{|\<ast\>|}>>w<rsub|i,c>
    log L<rsub|c><around*|(|X<rsub|i>\|\<theta\><rsub|c>|)>,>>>>
  </eqnarray*>

  where the first term does not contain the parameters, and the second term
  turns out to be the maximization objective in our M-step.\ 

  <subsection|The E-step and a marginal characterization of
  <math|q<around*|(|\<cdot\>|)>>>

  We show that the E-step, which solves for <math|<around*|{|w<rsub|i,c>|}>>,
  attempts to reduce <math|KL<around*|(|q<around*|\|||\|>P<around*|(|<around*|{|Z<rsub|i>|}>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>|)>>
  so that <math|q<around*|(|\<cdot\>|)>> better approximates the posterior of
  the latent variables. <math|q<around*|(|\<cdot\>|)>> is rewritten as

  <\equation*>
    q<around*|(|<around*|{|Z<rsub|i>|}>|)>=w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>\|Z<rsub|i>|)>\<approx\>w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>,
  </equation*>

  where an approximation is made by dropping the dependence on
  <math|Z<rsub|i>> in the conditional <math|q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>\|Z<rsub|i>|)>>.
  And <math|w<rsub|i><around*|(|Z<rsub|i>|)>=<big|prod><rsub|c>w<rsub|i,c><rsup|\<b-1\><around*|(|Z<rsub|i>=c|)>>>
  is a shorthand for the marginal of <math|Z<rsub|i>> under
  <math|q<around*|(|\<cdummy\>|)>>, subject to
  <math|<big|sum><rsub|c\<in\><around*|{|1,\<ldots\>,K|}>\<cup\><around*|{|\<ast\>|}>>w<rsub|i,c>=1>.

  Now we fix <math|i> and treat KL as a function of
  <math|<around*|{|w<rsub|i,c>|}>>,\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|KL<around*|(|q<around*|\|||\|>P<around*|(|<around*|{|Z<rsub|i>|}>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>|)>>|<cell|=>|<cell|<big|sum><rsub|Z<rsub|i>><big|sum><rsub|<around*|{|Z<rsub|-i>|}>>w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>
    log<frac|w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>|P<around*|(|<around*|{|Z<rsub|i>|}>,<around*|{|X<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>>+const,>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|Z<rsub|i>><big|sum><rsub|<around*|{|Z<rsub|-i>|}>>w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>
    log<frac|w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>|\<pi\><around*|(|Z<rsub|i>\|<around*|{|Z<rsub|-i>|}>|)>\<pi\><rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>L<rsub|Z<rsub|i>><around*|(|X<rsub|i>|)><big|prod><rsub|j\<neq\>i>L<rsub|Z<rsub|j>><around*|(|X<rsub|j>|)>>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|<around*|{|Z<rsub|i>|}>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log w<rsub|i><around*|(|Z<rsub|i>|)>-<big|sum><rsub|Z<rsub|i>><big|sum><rsub|<around*|{|Z<rsub|-i>|}>>w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>
    log <around*|(|<big|sum><rsub|k=1><rsup|K><big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=k|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=\<ast\>|)>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|)>-<big|sum><rsub|Z<rsub|i>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log <big|prod><rsub|c>L<rsub|c><around*|(|X<rsub|i>|)><rsup|<with|font-series|bold|1><around*|(|Z<rsub|i>=c|)>>+const>>>>
  </eqnarray*>

  where the const is due to the terms unrelated to
  <math|w<rsub|i><around*|(|Z<rsub|i>|)>>. The prior
  <math|\<pi\><around*|(|\<cdot\>|)>> is also decomposed into the conditional
  <math|\<pi\><around*|(|Z<rsub|i>\|<around*|{|Z<rsub|-i>|}>|)>> and the
  marginal <math|\<pi\><rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>>,
  where the conditional part can be expressed with the prediction rule
  (<reference|eqs:pik-prediction>). We make another approximation as follows
  <marginal-note|normal|c|This makes me uncomfortable because the
  <math|\<approx\>> below is really <math|\<leqslant\>>, which when inversed
  sign is a lower bound on KL. It certainly makes more sense to minimize an
  upper bound!>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|Z<rsub|i>><big|sum><rsub|<around*|{|Z<rsub|-i>|}>>w<rsub|i><around*|(|Z<rsub|i>|)>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)>
    log <around*|(|<big|sum><rsub|k=1><rsup|K><big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=k|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=\<ast\>|)>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|)>>>|<row|<cell|>|<cell|\<approx\>>|<cell|<big|sum><rsub|Z<rsub|i>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log <around*|(|<big|sum><rsub|<around*|{|Z<rsub|-i>|}>>q<rprime|'><around*|(|<around*|{|Z<rsub|-i>|}>|)><big|sum><rsub|k=1><rsup|K><big|sum><rsub|j\<neq\>i><with|font-series|bold|1><around*|(|Z<rsub|j>=k|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>\<b-1\><around*|(|Z<rsub|j>=\<ast\>|)>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|Z<rsub|i>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log <around*|(|<big|sum><rsub|k=1><rsup|K><big|sum><rsub|j\<neq\>i>\<bbb-E\><rsub|q<rprime|'>>\<b-1\><around*|(|Z<rsub|j>=k|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>\<bbb-E\><rsub|q<rprime|'>>\<b-1\><around*|(|Z<rsub|j>=\<ast\>|)>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|Z<rsub|i>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log <around*|[|<big|sum><rsub|k=1><rsup|K><around*|(|<big|sum><rsub|j\<neq\>i>w<rsub|j,k>|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>w<rsub|j,\<ast\>>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|]>.>>>>
  </eqnarray*>

  Now, we have\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|KL<around*|(|q<around*|\|||\|>P<around*|(|<around*|{|Z<rsub|i>|}>\|<around*|{|X<rsub|i>|}>,<around*|{|\<theta\><rsub|c>|}>|)>|)>>>|<row|<cell|>|<cell|\<approx\>>|<cell|<big|sum><rsub|Z<rsub|i>>w<rsub|i><around*|(|Z<rsub|i>|)>
    log <frac|w<rsub|i><around*|(|Z<rsub|i>|)>|<around*|[|<big|sum><rsub|k=1><rsup|K><around*|(|<big|sum><rsub|j\<neq\>i>w<rsub|j,k>|)>\<delta\><rsub|k><around*|(|Z<rsub|i>|)>+<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>w<rsub|j,\<ast\>>|)>\<delta\><rsub|\<ast\>><around*|(|Z<rsub|i>|)>|]><big|prod><rsub|c>L<rsub|c><around*|(|X<rsub|i>|)><rsup|\<b-1\><around*|(|Z<rsub|i>=c|)>>>+const,>>>>
  </eqnarray*>

  which takes the form of another KL divergence. Clearly, this is minimized
  by matching the two distributions on <math|Z<rsub|i>>, which gives the
  M-step equations\ 

  <\equation*>
    w<rsub|i,c>\<propto\><choice|<tformat|<table|<row|<cell|<around*|(|\<alpha\>+<big|sum><rsub|j\<neq\>i>w<rsub|j,\<ast\>>|)><big|int>P<around*|(|X<rsub|i>\|\<theta\><rsup|\<ast\>>|)>G<rsub|0><around*|(|\<theta\><rsup|\<ast\>>|)>\<mathd\>\<theta\><rsup|\<ast\>><space|1em><around*|(|c=\<ast\>|)>>>|<row|<cell|<around*|(|<big|sum><rsub|j\<neq\>i>w<rsub|j,c>|)>P<around*|(|X<rsub|i><mid|\|>\<theta\><rsub|c>|)><space|1em><around*|(|1\<leqslant\>c\<leqslant\>K|)>>>>>>,
  </equation*>

  subject to <math|<big|sum><rsub|c>w<rsub|i,c>=1>. Previously in the
  description of the algorithm, we intuitively interpreted
  <math|<around*|(|<big|sum><rsub|j\<neq\>i>w<rsub|j,c>|)>> as an estimate of
  the size of cluster <math|c>.

  <subsection|Justifying the Bayesian posterior predictive variant>

  In the Bayesian posterior predictive variant of the algorithm, in the
  M-step, we instead update the estimate of parameters in terms of its
  posterior distribution. Supposing the posterior for
  <math|\<theta\><rsub|c>> is <math|H<around*|(|\<theta\><rsub|c>\|\<psi\><rsub|c>|)>>,
  it is updated to\ 

  <\equation*>
    H<around*|(|\<theta\><rsub|c>\|\<psi\><rsub|c>|)>=G<rsub|0><around*|(|\<theta\><rsub|c>|)><big|prod><rsub|i=1><rsup|n>P<around*|(|X<rsub|i><mid|\|>\<theta\><rsub|c>|)><rsup|w<rsub|i,c>>.
  </equation*>

  Why is <math|<big|prod><rsub|i=1><rsup|n>P<around*|(|X<rsub|i>\|\<theta\><rsub|c>|)><rsup|w<rsub|i,c>>>
  the right likelihood to be used? In (<reference|eqs:EM>), when
  <math|q<around*|(|\<cdot\>|)>> closely approximates the true posterior, we
  have\ 

  <\equation*>
    log P<around*|(|<around*|{|X<rsub|i>|}>\|<around*|{|\<theta\><rsub|c>|}>|)>\<approx\>\<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>+H<around*|(|q|)>=\<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>+const,
  </equation*>

  where\ 

  <\equation*>
    \<cal-L\><around*|(|<around*|{|\<theta\><rsub|c>|}>|)>=<big|sum><rsub|<around*|{|Z<rsub|i>|}>>q<around*|(|<around*|{|Z<rsub|i>|}>|)>
    log \<pi\><around*|(|<around*|{|Z<rsub|i>|}>|)>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|c\<in\><around*|[|K|]>\<cup\><around*|{|\<ast\>|}>>w<rsub|i,c>
    log L<rsub|c><around*|(|X<rsub|i>\|\<theta\><rsub|c>|)>,
  </equation*>

  which relates to <math|\<theta\><rsub|c>> via
  <math|<big|sum><rsub|i=1><rsup|n>w<rsub|i,c> log
  L<rsub|c><around*|(|X<rsub|i>\|\<theta\><rsub|c>|)>>, which becomes
  <math|<big|prod><rsub|i=1><rsup|n>P<around*|(|X<rsub|i>\|\<theta\><rsub|c>|)><rsup|w<rsub|i,c>>>
  when exponentiated.

  \;
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
    <associate|page-type|letter>
    <associate|page-width|auto>
    <associate|page-width-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|3|2>>
    <associate|auto-5|<tuple|3.1|3>>
    <associate|auto-6|<tuple|3.2|4>>
    <associate|auto-7|<tuple|3.3|5>>
    <associate|auto-8|<tuple|4|5>>
    <associate|eqs:EM|<tuple|1|3>>
    <associate|eqs:pik-prediction|<tuple|1|2>>
    <associate|fig:pik|<tuple|1|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|An example where a draw from CRP is collapsed into a
      <with|mode|<quote|math>|\<pi\><rsub|2>> draw. |<pageref|auto-2>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Collapsing
      CRP into finite-dimensional <with|mode|<quote|math>|\<pi\><rsub|K>>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Theoretical
      Overview> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Theory
      for the EM procedure> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>The M-Step
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>The E-step and a marginal
      characterization of <with|mode|<quote|math>|q<around*|(|\<cdot\>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Justifying the Bayesian
      posterior predictive variant <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>
    </associate>
  </collection>
</auxiliary>