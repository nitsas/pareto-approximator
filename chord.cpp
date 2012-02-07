/* chord.cpp */


// Won't #include "chord.h". In fact "chord.h" will include "chord.cpp" 
// because it describes function templates (which doesn't allow us to 
// split declaration from implementation.


using std::min;


namespace pareto_approximator {


template <class S> 
list< PointAndSolution<S> > 
chordAlgorithm(function<PointAndSolution<S> (double, double)> comb, double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

  // find the westmost (best on x objective) and the southmost (best on y
  // objective) solutions (and the corresponding points in objective space)
  PointAndSolution<S> west = comb(1.0, 0.0);
  PointAndSolution<S> south = comb(0.0, 1.0);

  // find that point in objective space which corresponds to the best possible 
  // solution we could expect (we'll pass it to doChord as a sort of lower 
  // limit on the search space)
  Point tip(min(west.point.x, south.point.x), min(west.point.y, south.point.y));

  // let doChord do all the work (it's recursive)
  return doChord<S>(comb, west, south, tip, eps);
}


template <class S> 
list< PointAndSolution<S> > 
doChord(function<PointAndSolution<S> (double, double)> comb, 
        const PointAndSolution<S>& west, const PointAndSolution<S>& south, 
        const Point& tip, double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

  // check if the best possible point is approximately dominated by the 
  // points we have so far
  Line2D ws(west.point, south.point);
  if (ws.ratioDistance(tip) <= eps) {
    list< PointAndSolution<S> > resultList;
    resultList.push_back(west);
    resultList.push_back(south);
    return resultList;
  }
  // else

  // check if ws is a vertical line (slope is infinite) and call comb 
  // using its slope
  PointAndSolution<S> southwest;
  if (!ws.isVertical())
    southwest = comb(-ws.m(), 1.0);
  else
    southwest = comb(1.0, 0.0);

  // check if the point we just found is approximately dominated by the 
  // points we have so far
  if (ws.ratioDistance(southwest.point) <= eps) {
    list< PointAndSolution<S> > resultList;
    resultList.push_back(west);
    resultList.push_back(south);
    return resultList;
  }
  // else

  // split the problem into two subproblems, the west one and the south one
  Line2D parallel = ws.parallelThrough(southwest.point);
  Line2D wt(west.point, tip);
  Line2D ts(tip, south.point);
  Point westTip  = parallel.intersection(wt);
  Point southTip = parallel.intersection(ts);
  list< PointAndSolution<S> > westList, southList;
  // call comb on the two subproblems
  westList  = doChord<S>(comb, west,  southwest, westTip,  eps);
  southList = doChord<S>(comb, southwest, south, southTip, eps);
  // remove southList's first element (it's the same as westList's last one)
  // and merge the lists
  southList.pop_front();
  westList.merge(southList);

  return westList;
}


}  // namespace pareto_approximator


