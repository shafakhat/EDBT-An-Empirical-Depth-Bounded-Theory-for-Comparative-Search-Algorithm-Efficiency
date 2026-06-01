
/**
 * Represents a brief description of EDBT: An Accuracy-Aware Extension of Behavioural Tree Search.
 * @Shafakhatullah Khan Mohammed (shafakhath91@gmail.com)
 * @version 1.0
 * @since 2025-10-18
 */

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

public class UnifiedTreeExperimentL extends JPanel {

    static final int    NUM_TRIALS    = 30;
    static final int[]  DEPTHS        = {4, 6, 8, 10, 12};
    static final int    BRANCHING     = 4;
    static final double OBSTACLE_RATE = 0.15;
    static final int    AB_MAX_DEPTH  = 9;
    static final double EPSILON       = 0.05;

    // ============================================================
    // RESULT RECORD
    // ============================================================
    static class AlgResult {
        String  name;
        int     depth;
        double  meanNodes, stdNodes, meanUseful, meanBEI, solutionRate;
        double  meanDeltaD, meanDeltaS, meanDeltaC, meanABEI, meanQstar;
        double  certDeltaS, acpProduct, gammaMinus, gammaPlus, gammaMu;
        boolean acpSatisfied, pcctSatisfied;

        AlgResult(String n, int d, double mn, double sn,
                  double mu, double mb, double sr,
                  double dd, double ds, double dc, double abei,
                  double qs, double cds, double acp,
                  double gm, double gp, double gmu) {
            name = n; depth = d;
            meanNodes = mn; stdNodes = sn; meanUseful = mu;
            meanBEI = mb; solutionRate = sr;
            meanDeltaD = dd; meanDeltaS = ds; meanDeltaC = dc;
            meanABEI = abei; meanQstar = qs; certDeltaS = cds;
            acpProduct = acp; gammaMinus = gm; gammaPlus = gp; gammaMu = gmu;
            acpSatisfied  = Double.isFinite(acp) && Double.isFinite(gm) && (acp >= gm - 1e-9);
            pcctSatisfied = Double.isFinite(ds)  && Double.isFinite(cds) && (ds  >= cds - 1e-9);
        }
    }

    // ============================================================
    // EDBT COMPONENT 1: Complexity Profile
    // ============================================================
    static double clampBeta(double beta) {
        if (!Double.isFinite(beta) || beta < 0.01) return 0.01;
        if (beta > 1.0) return 1.0;
        return beta;
    }

    static double[] complexityProfile(double betaRaw, int d, int b, double rho) {
        double beta   = clampBeta(betaRaw);
        double bBeta  = Math.pow(b, beta);
        double denom  = bBeta - 1.0;
        if (Math.abs(denom) < 1e-10) denom = 1e-10;
        double alpha  = bBeta / denom;
        double bBetaD = Math.pow(b, beta * d);
        double gammaPlus  = alpha * bBetaD;
        double kmDenom    = 1.0 - Math.pow(b, -beta / 2.0);
        if (Math.abs(kmDenom) < 1e-10) kmDenom = 1e-10;
        double gammaMinus = Math.pow(b, beta * d / 2.0) / kmDenom;
        double gammaMu    = alpha * Math.pow(1.0 - rho, d) * bBetaD;
        gammaPlus  = Double.isFinite(gammaPlus)  && gammaPlus  > 0 ? gammaPlus  : 1.0;
        gammaMinus = Double.isFinite(gammaMinus) && gammaMinus > 0 ? gammaMinus : 1.0;
        gammaMu    = Double.isFinite(gammaMu)    && gammaMu    > 0 ? gammaMu    : 1.0;
        return new double[]{gammaPlus, gammaMinus, gammaMu};
    }

    // ============================================================
    // EDBT COMPONENT 2: Accuracy Function
    // ============================================================
    static double[] accuracyMeasures(double[] psiAll, double[] psiKept) {
        if (psiAll == null || psiAll.length == 0)
            return new double[]{1.0, 0.0};
        if (psiKept == null || psiKept.length == 0)
            return new double[]{-1.0, -1.0};
        if (psiKept.length == psiAll.length)
            return new double[]{1.0, 0.0};
        double minAll = Double.MAX_VALUE;
        for (double v : psiAll) if (v < minAll) minAll = v;
        double minKept = Double.MAX_VALUE;
        for (double v : psiKept) if (v < minKept) minKept = v;
        double deltaS = (minAll > 1e-12) ? Math.min(1.0, minKept / minAll) : 1.0;
        double deltaC = 1.0 - (double) psiKept.length / (double) psiAll.length;
        return new double[]{deltaS, deltaC};
    }

    // ============================================================
    // EDBT COMPONENT 3: Convergence Operator
    // ============================================================
    static double[] convergenceOperator(double betaRaw, double deltaD,
                                        double eps, int b) {
        double beta  = clampBeta(betaRaw);
        double bBeta = Math.pow(b, beta);
        double denom = bBeta - 1.0;
        if (Math.abs(denom) < 1e-10) denom = 1e-10;
        double alpha = bBeta / denom;
        double dStar;
        if (eps <= 0 || eps >= 1 || alpha <= 0 || !Double.isFinite(alpha)) {
            dStar = Double.MAX_VALUE;
        } else {
            double arg = (1.0 - eps) / (eps * alpha);
            if (arg <= 0 || !Double.isFinite(arg)) {
                dStar = Double.MAX_VALUE;
            } else if (arg <= 1.0) {
                dStar = 1.0;
            } else {
                dStar = Math.ceil((1.0 / beta) * (Math.log(arg) / Math.log(b)));
                if (dStar < 1) dStar = 1;
            }
        }
        double qStar = deltaD * (1.0 - eps);
        return new double[]{dStar, qStar};
    }

    // ============================================================
    // EDBT: Path-level deltaD
    // ============================================================
    static double computePathDeltaD(Grid g, int[] par, int start, int goal,
                                    double wt, boolean isAdmissible) {
        if (par == null || (par[goal] == -1 && goal != start)) return 1.0;
        List<Integer> pathNodes = new ArrayList<>();
        int cur = goal;
        Set<Integer> seen = new HashSet<>();
        while (cur != start) {
            if (seen.contains(cur)) return 1.0;
            seen.add(cur);
            pathNodes.add(cur);
            if (par[cur] < 0) return 1.0;
            cur = par[cur];
        }
        pathNodes.add(start);
        Collections.reverse(pathNodes);
        if (pathNodes.size() < 2) return 1.0;
        double deltaD  = 1.0;
        double[] gCost = new double[g.total()];
        Arrays.fill(gCost, Double.MAX_VALUE);
        gCost[start] = 0;
        for (int i = 0; i < pathNodes.size() - 1; i++) {
            int u = pathNodes.get(i);
            int v = pathNodes.get(i + 1);
            int edgeCost = 1;
            for (int[] nb : g.nbrs(u))
                if (g.id(nb[0], nb[1]) == v) { edgeCost = nb[2]; break; }
            gCost[v] = gCost[u] + edgeCost;
        }
        for (int i = 0; i < pathNodes.size() - 1; i++) {
            int u = pathNodes.get(i);
            List<int[]> nbrs = g.nbrs(u);
            if (nbrs.isEmpty()) continue;
            double[] psiAll = new double[nbrs.size()];
            for (int j = 0; j < nbrs.size(); j++) {
                int nid = g.id(nbrs.get(j)[0], nbrs.get(j)[1]);
                psiAll[j] = gCost[u] + nbrs.get(j)[2] + wt * g.h(nid);
            }
            int    nextNode = pathNodes.get(i + 1);
            double gNext    = gCost[nextNode];
            double psiNext  = gNext + wt * g.h(nextNode);
            double[] acc    = accuracyMeasures(psiAll, new double[]{psiNext});
            if (acc[0] >= 0) deltaD *= acc[0];
        }
        return deltaD;
    }

    // ============================================================
    // PCCT Certificate
    // ============================================================
    static double pcctCertificate(String algName) {
        if (algName.startsWith("WA*(1.5)")) return 1.0 / 1.5;
        if (algName.startsWith("WA*(2.0)")) return 1.0 / 2.0;
        if (algName.startsWith("WA*(5.0)")) return 1.0 / 5.0;
        if (algName.equals("GBFS"))         return 0.0;
        return 1.0;
    }

    // ============================================================
    // ABEI and BEI
    // ============================================================
    static double computeABEI(double N, double U, double q, double deltaD, int b) {
        if (N < 1 || U < 1) return 0.0;
        double logNorm = Math.log(N + 1) / Math.log(b);
        if (logNorm < 1e-12) return 0.0;
        return (U * q * deltaD) / (N * logNorm);
    }

    static double computeBEI(double N, double U, double q, int b) {
        if (N < 1 || U < 1) return 0.0;
        double logNorm = Math.log(N + 1) / Math.log(b);
        if (logNorm < 1e-12) return 0.0;
        return (U * q) / (N * logNorm);
    }

    // ============================================================
    // GRID
    // ============================================================
    static class Grid {
        final int R, C;
        final int[][] w;
        final boolean[][] blocked;
        final int sr, sc, gr, gc;
        static final int[] DR = {-1, 1,  0, 0};
        static final int[] DC = { 0, 0, -1, 1};

        Grid(int size, long seed) {
            R = C = size;
            w       = new int[R][C];
            blocked = new boolean[R][C];
            sr = sc = 0; gr = R - 1; gc = C - 1;
            Random rng = new Random(seed);
            for (int r = 0; r < R; r++)
                for (int c = 0; c < C; c++)
                    w[r][c] = 1 + rng.nextInt(10);
            for (int r = 0; r < R; r++)
                for (int c = 0; c < C; c++) {
                    if ((r == sr && c == sc) || (r == gr && c == gc)) continue;
                    if (rng.nextDouble() < OBSTACLE_RATE) blocked[r][c] = true;
                }
            for (int c = 0; c < C; c++) blocked[0][c]     = false;
            for (int r = 0; r < R; r++) blocked[r][C - 1] = false;
        }

        boolean valid(int r, int c) { return r>=0&&r<R&&c>=0&&c<C&&!blocked[r][c]; }
        int id(int r, int c)  { return r * C + c; }
        int row(int id)       { return id / C; }
        int col(int id)       { return id % C; }
        int total()           { return R * C; }
        int h(int r, int c)   { return Math.abs(r - gr) + Math.abs(c - gc); }
        int h(int id)         { return h(row(id), col(id)); }

        List<int[]> nbrs(int r, int c) {
            List<int[]> out = new ArrayList<>(4);
            for (int d = 0; d < 4; d++) {
                int nr = r + DR[d], nc = c + DC[d];
                if (valid(nr, nc)) out.add(new int[]{nr, nc, w[nr][nc]});
            }
            return out;
        }
        List<int[]> nbrs(int id) { return nbrs(row(id), col(id)); }
    }

    // ============================================================
    // PATH LENGTH
    // ============================================================
    static int pathLengthArr(int[] parent, int start, int goal, int maxLen) {
        if (goal == start) return 1;
        if (parent == null || parent[goal] == -1) return 0;
        int len = 1, cur = goal;
        Set<Integer> vis = new HashSet<>();
        while (cur != start) {
            if (vis.contains(cur) || len > maxLen) return 0;
            vis.add(cur); int p = parent[cur];
            if (p < 0) return 0;
            cur = p; len++;
        }
        return len;
    }

    // ============================================================
    // EDBT METRICS
    // ============================================================
    static class EDBTMetrics {
        long    nodes;
        int     useful;
        double  deltaD;
        double  deltaSMean;
        double  deltaCMean;
        int[]   par;
        boolean found;
    }

    // ============================================================
    // ALGORITHM 1: BFS
    // ============================================================
    static EDBTMetrics bfsExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        int[]     par = new int[T];     Arrays.fill(par, -1);
        boolean[] vis = new boolean[T];
        Queue<Integer> q = new LinkedList<>();
        q.add(start); vis[start] = true;
        long expanded = 0; boolean found = false;
        while (!q.isEmpty()) {
            int u = q.poll(); expanded++;
            if (u == goal) { found = true; break; }
            for (int[] nb : g.nbrs(u)) {
                int nid = g.id(nb[0], nb[1]);
                if (!vis[nid]) { vis[nid] = true; par[nid] = u; q.add(nid); }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found;
        m.deltaD = 1.0; m.deltaSMean = 1.0; m.deltaCMean = 0.0; m.par = par;
        return m;
    }

    // ============================================================
    // ALGORITHM 2: DIJKSTRA
    // ============================================================
    static EDBTMetrics dijkstraExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        double[] dist = new double[T]; Arrays.fill(dist, Double.MAX_VALUE);
        int[]    par  = new int[T];    Arrays.fill(par, -1);
        dist[start] = 0;
        PriorityQueue<double[]> pq = new PriorityQueue<>(Comparator.comparingDouble(x -> x[0]));
        pq.add(new double[]{0, start});
        boolean[] vis = new boolean[T];
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!pq.isEmpty()) {
            double[] cur = pq.poll(); int u = (int) cur[1];
            if (vis[u]) continue;
            vis[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = dist[u] + nbrs.get(i)[2];
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]); double nd = dist[u] + nb[2];
                    if (!vis[nid] && nd < dist[nid]) kl.add(nd);
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : g.nbrs(u)) {
                int nid = g.id(nb[0], nb[1]); double nd = dist[u] + nb[2];
                if (nd < dist[nid]) { dist[nid] = nd; par[nid] = u; pq.add(new double[]{nd, nid}); }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, 1.0, true) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 3/4: A* / WA*
    // ============================================================
    static EDBTMetrics astarExpand(Grid g, double wt) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        double[] gCost = new double[T]; Arrays.fill(gCost, Double.MAX_VALUE);
        int[]    par   = new int[T];    Arrays.fill(par, -1);
        gCost[start] = 0;
        PriorityQueue<double[]> open = new PriorityQueue<>(Comparator.comparingDouble(x -> x[0]));
        open.add(new double[]{wt * g.h(start), start});
        boolean[] closed = new boolean[T];
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!open.isEmpty()) {
            double[] cur = open.poll(); int u = (int) cur[1];
            if (closed[u]) continue;
            closed[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++) {
                    int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                    psiAll[i] = gCost[u] + nbrs.get(i)[2] + wt * g.h(nid);
                }
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]); double ng = gCost[u] + nb[2];
                    if (!closed[nid] && ng < gCost[nid]) kl.add(ng + wt * g.h(nid));
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (closed[nid]) continue;
                double ng = gCost[u] + nb[2];
                if (ng < gCost[nid]) {
                    gCost[nid] = ng; par[nid] = u;
                    open.add(new double[]{ng + wt * g.h(nid), nid});
                }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, wt, wt <= 1.0) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 5: GBFS
    // ============================================================
    static EDBTMetrics gbfsExpand(Grid g) {
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        int T = g.total();
        int[]     par    = new int[T];     Arrays.fill(par, -1);
        boolean[] vis    = new boolean[T];
        boolean[] inOpen = new boolean[T];
        PriorityQueue<int[]> open = new PriorityQueue<>(Comparator.comparingInt(x -> x[0]));
        open.add(new int[]{g.h(start), start}); inOpen[start] = true;
        long expanded = 0; boolean found = false;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (!open.isEmpty()) {
            int[] cur = open.poll(); int u = cur[1];
            if (vis[u]) continue;
            vis[u] = true; expanded++;
            if (u == goal) { found = true; break; }
            List<int[]> nbrs = g.nbrs(u);
            if (!nbrs.isEmpty()) {
                double[] psiAll = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = g.h(g.id(nbrs.get(i)[0], nbrs.get(i)[1]));
                List<Double> kl = new ArrayList<>();
                for (int[] nb : nbrs) {
                    int nid = g.id(nb[0], nb[1]);
                    if (!vis[nid]) kl.add((double) g.h(nid));
                }
                if (!kl.isEmpty()) {
                    double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                    if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                }
            }
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (!vis[nid] && !inOpen[nid]) {
                    if (par[nid] == -1) par[nid] = u;
                    open.add(new int[]{g.h(nid), nid});
                    inOpen[nid] = true;
                }
            }
        }
        int useful = found ? pathLengthArr(par, start, goal, T) : 0;
        double deltaD = found ? computePathDeltaD(g, par, start, goal, 0.0, false) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = found; m.par = par;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 6: BIDIRECTIONAL DIJKSTRA
    // ============================================================
    static EDBTMetrics biDijkstraExpand(Grid g) {
        int start = g.id(g.sr,g.sc), goal = g.id(g.gr,g.gc), T = g.total();
        double[] dF = new double[T], dB = new double[T];
        Arrays.fill(dF, Double.MAX_VALUE); Arrays.fill(dB, Double.MAX_VALUE);
        dF[start] = 0; dB[goal] = 0;
        boolean[] vF = new boolean[T], vB = new boolean[T];
        int[] parF = new int[T], parB = new int[T];
        Arrays.fill(parF, -1); Arrays.fill(parB, -1);
        PriorityQueue<double[]> pqF = new PriorityQueue<>(Comparator.comparingDouble(x->x[0]));
        PriorityQueue<double[]> pqB = new PriorityQueue<>(Comparator.comparingDouble(x->x[0]));
        pqF.add(new double[]{0, start}); pqB.add(new double[]{0, goal});
        long expanded = 0; double mu = Double.MAX_VALUE; int meet = -1;
        double sumDS = 0, sumDC = 0; int steps = 0;
        while (true) {
            if (pqF.isEmpty() && pqB.isEmpty()) break;
            double topF = pqF.isEmpty() ? Double.MAX_VALUE : pqF.peek()[0];
            double topB = pqB.isEmpty() ? Double.MAX_VALUE : pqB.peek()[0];
            if (topF + topB >= mu) break;
            if (!pqF.isEmpty() && topF <= topB) {
                double[] cur = pqF.poll(); int u = (int) cur[1];
                if (vF[u]) continue;
                vF[u] = true; expanded++;
                List<int[]> nbrs = g.nbrs(u);
                if (!nbrs.isEmpty()) {
                    double[] psiAll = new double[nbrs.size()];
                    for (int i = 0; i < nbrs.size(); i++) psiAll[i] = dF[u] + nbrs.get(i)[2];
                    List<Double> kl = new ArrayList<>();
                    for (int[] nb : nbrs) {
                        int nid = g.id(nb[0], nb[1]); double nd = dF[u] + nb[2];
                        if (nd < dF[nid]) kl.add(nd);
                    }
                    if (!kl.isEmpty()) {
                        double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                        if (acc[0] >= 0) { sumDS += acc[0]; sumDC += acc[1]; steps++; }
                    }
                }
                for (int[] nb : g.nbrs(u)) {
                    int nid = g.id(nb[0], nb[1]); double nd = dF[u] + nb[2];
                    if (nd < dF[nid]) { dF[nid] = nd; parF[nid] = u; pqF.add(new double[]{nd, nid}); }
                    if (vB[nid] && dF[u]+nb[2]+dB[nid] < mu) { mu = dF[u]+nb[2]+dB[nid]; meet = nid; }
                }
            } else if (!pqB.isEmpty()) {
                double[] cur = pqB.poll(); int u = (int) cur[1];
                if (vB[u]) continue;
                vB[u] = true; expanded++;
                for (int[] nb : g.nbrs(u)) {
                    int nid = g.id(nb[0], nb[1]); double nd = dB[u] + nb[2];
                    if (nd < dB[nid]) { dB[nid] = nd; parB[nid] = u; pqB.add(new double[]{nd, nid}); }
                    if (vF[nid] && dB[u]+nb[2]+dF[nid] < mu) { mu = dB[u]+nb[2]+dF[nid]; meet = nid; }
                }
            } else break;
        }
        int useful = 0;
        if (meet >= 0) {
            useful  = pathLengthArr(parF, start, meet, T);
            useful += pathLengthArr(parB, goal,  meet, T) - 1;
        }
        double deltaD = (meet >= 0 && useful > 0)
                ? computePathDeltaD(g, parF, start, meet, 1.0, true) : 1.0;
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = expanded; m.useful = useful; m.found = (meet >= 0); m.par = parF;
        m.deltaD = deltaD;
        m.deltaSMean = (steps > 0) ? sumDS / steps : 1.0;
        m.deltaCMean = (steps > 0) ? sumDC / steps : 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 7: RBFS
    // ============================================================
    static long   rbfsCount;
    static double rbfsSumDS;
    static int    rbfsSteps;

    // ---- RBFS result container to avoid negative-return ambiguity ----
    static boolean rbfsFound;
    static double  rbfsBestCost;

    /**
     * Returns: f-value to propagate up (always >= 0 and finite, or MAX_VALUE).
     * Goal detection is signaled via rbfsFound flag (avoids negative/MAX confusion).
     */
    static double rbfsRec(Grid g, int cur, double gCur, double fCur,
                           double fLimit, int goal,
                           Set<Integer> onPath,
                           List<int[]> pathEdges,
                           Map<Integer, Double> bestG) {  

        // Hard cap checked immediately at entry
        if (rbfsCount >= 2_000_000) {
            rbfsFound = false;
            return Double.MAX_VALUE;
        }
        rbfsCount++;

        // Goal check
        if (cur == goal) {
            rbfsFound    = true;
            rbfsBestCost = gCur;
            return 0.0; // signals "found" — caller checks rbfsFound
        }

        // Generate successors with bestG pruning
        List<double[]> succs = new ArrayList<>();
        for (int[] nb : g.nbrs(cur)) {
            int    nid      = g.id(nb[0], nb[1]);
            double edgeCost = nb[2];
            double ng       = gCur + edgeCost;

            // Cycle guard (onPath = current DFS path)
            if (onPath.contains(Integer.valueOf(nid))) continue;

            // BestG pruning — only prune if not strictly better
            Double knownG = bestG.get(nid);
            if (knownG != null && ng >= knownG) continue;

            // Record best g (will be rolled back on backtrack if needed)
            bestG.put(nid, ng);

            double nf = Math.max(ng + g.h(nid), fCur);
            succs.add(new double[]{nid, edgeCost, nf});
        }

        // Dead end
        if (succs.isEmpty()) return Double.MAX_VALUE;

        // EDBT accuracy measurement
        double[] psiAll = succs.stream().mapToDouble(s -> s[2]).toArray();
        succs.sort(Comparator.comparingDouble(s -> s[2]));
        double[] psiKept = new double[]{succs.get(0)[2]};
        double[] acc     = accuracyMeasures(psiAll, psiKept);
        if (acc[0] >= 0) { rbfsSumDS += acc[0]; rbfsSteps++; }

        // Main RBFS loop
        // Guard against infinite loop when fLimit = MAX_VALUE
        int maxIter = succs.size() * 2_000_000; // generous but bounded
        int iter    = 0;

        while (iter++ < maxIter) {

            // Re-sort to get current best
            succs.sort(Comparator.comparingDouble(s -> s[2]));
            double[] best = succs.get(0);

            // Explicit MAX_VALUE guard
            if (best[2] >= Double.MAX_VALUE) return Double.MAX_VALUE;

            // Standard f-limit check
            if (best[2] > fLimit) return best[2];

            double altF    = succs.size() > 1 ? succs.get(1)[2] : Double.MAX_VALUE;
            int    bestNid = (int) best[0];
            double bestEC  = best[1];

            // Expand best successor
            onPath.add(Integer.valueOf(bestNid));
            pathEdges.add(new int[]{cur, bestNid, (int) bestEC});

            double newLimit = (altF >= Double.MAX_VALUE) ? fLimit : Math.min(fLimit, altF);
            double res = rbfsRec(g, bestNid, gCur + bestEC, best[2],
                                  newLimit, goal, onPath, pathEdges, bestG);

            // Check if goal was found (via flag — avoids return-value ambiguity)
            if (rbfsFound) return 0.0; // propagate found upward

            // Backtrack
            pathEdges.remove(pathEdges.size() - 1);
            onPath.remove(Integer.valueOf(bestNid));

            // Update f-value with returned bound
            // If res is MAX_VALUE, mark this successor as exhausted
            best[2] = (res >= Double.MAX_VALUE) ? Double.MAX_VALUE : res;
        }

        // Exhausted iterations safely
        return Double.MAX_VALUE;
    }

    static EDBTMetrics rbfsExpand(Grid g) {
        rbfsCount = 0; rbfsSumDS = 0.0; rbfsSteps = 0;
        rbfsFound = false; rbfsBestCost = Double.MAX_VALUE;

        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);

        // BestG is local per search — no cross-trial contamination
        Map<Integer, Double> bestG = new HashMap<>();
        bestG.put(start, 0.0);

        Set<Integer> onPath    = new HashSet<>();
        List<int[]>  pathEdges = new ArrayList<>();
        onPath.add(Integer.valueOf(start));

        // Run RBFS with increasing f-limits (like IDA* outer loop)
        // This ensures termination even if single call hits fLimit
        double fLimit = g.h(start);
        boolean done  = false;

        for (int iter = 0; iter < 100_000 && !done; iter++) {
            rbfsFound = false;
            pathEdges.clear();
            onPath.clear();
            onPath.add(Integer.valueOf(start));
            bestG.clear();
            bestG.put(start, 0.0);

            double res = rbfsRec(g, start, 0, fLimit, fLimit,
                                  goal, onPath, pathEdges, bestG);

            if (rbfsFound) {
                done = true;
            } else if (res >= Double.MAX_VALUE || rbfsCount >= 2_000_000) {
                break; // no solution or budget exhausted
            } else {
                fLimit = res; // increase threshold (IDA*-style)
            }
        }

        boolean found  = rbfsFound;
        int     useful = found ? pathEdges.size() + 1 : 0;

        // Compute path-level deltaD from pathEdges
        double deltaD = 1.0;
        if (found && !pathEdges.isEmpty()) {
            double gCur = 0;
            for (int[] edge : pathEdges) {
                int u = edge[0], v = edge[1], ec = edge[2];
                List<int[]> nbrs = g.nbrs(u);
                double[] psiAll  = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++) {
                    int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                    psiAll[i] = gCur + nbrs.get(i)[2] + g.h(nid);
                }
                double gNext    = gCur + ec;
                double[] psiKpt = new double[]{gNext + g.h(v)};
                double[] a      = accuracyMeasures(psiAll, psiKpt);
                if (a[0] >= 0) deltaD *= a[0];
                gCur = gNext;
            }
        }

        EDBTMetrics m = new EDBTMetrics();
        m.nodes      = rbfsCount;
        m.useful     = useful;
        m.found      = found;
        m.par        = null;
        m.deltaD     = deltaD;
        m.deltaSMean = (rbfsSteps > 0) ? rbfsSumDS / rbfsSteps : 1.0;
        m.deltaCMean = 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 8: IDA*
    // ============================================================
    static long   idaCount;
    static double idaSumDS;
    static int    idaStepsCount;

    static double idaSearch(Grid g, int cur, double gCur, double thresh,
                             int goal, Set<Integer> onPath, List<int[]> pathEdges) {
        double f = gCur + g.h(cur);
        if (f > thresh) return f;
        idaCount++;
        if (idaCount > 10_000_000) return Double.MAX_VALUE;
        if (cur == goal) return -gCur;

        List<int[]> nbrs = g.nbrs(cur);
        nbrs.sort(Comparator.comparingInt(nb -> (int)(gCur + nb[2]) + g.h(nb[0], nb[1])));

        if (!nbrs.isEmpty()) {
            double[] psiAll = new double[nbrs.size()];
            for (int i = 0; i < nbrs.size(); i++) {
                int nid = g.id(nbrs.get(i)[0], nbrs.get(i)[1]);
                psiAll[i] = gCur + nbrs.get(i)[2] + g.h(nid);
            }
            List<Double> kl = new ArrayList<>();
            for (int[] nb : nbrs) {
                int nid = g.id(nb[0], nb[1]);
                if (!onPath.contains(Integer.valueOf(nid)))
                    kl.add(gCur + nb[2] + (double) g.h(nid));
            }
            if (!kl.isEmpty()) {
                double[] acc = accuracyMeasures(psiAll, kl.stream().mapToDouble(x->x).toArray());
                if (acc[0] >= 0) { idaSumDS += acc[0]; idaStepsCount++; }
            }
        }

        double minExc = Double.MAX_VALUE;
        for (int[] nb : nbrs) {
            int nid = g.id(nb[0], nb[1]);
            if (onPath.contains(Integer.valueOf(nid))) continue;
            onPath.add(Integer.valueOf(nid));
            pathEdges.add(new int[]{cur, nid, nb[2]});
            double res = idaSearch(g, nid, gCur + nb[2], thresh, goal, onPath, pathEdges);
            if (res < 0) return res;
            pathEdges.remove(pathEdges.size() - 1);
            onPath.remove(Integer.valueOf(nid));
            if (res < minExc) minExc = res;
        }
        return minExc;
    }

    static EDBTMetrics idaExpand(Grid g) {
        idaCount = 0; idaSumDS = 0.0; idaStepsCount = 0;
        int start = g.id(g.sr, g.sc), goal = g.id(g.gr, g.gc);
        double thresh = g.h(start);
        Set<Integer> onPath    = new HashSet<>();
        List<int[]>  pathEdges = new ArrayList<>();
        onPath.add(Integer.valueOf(start));
        boolean found = false;
        for (int iter = 0; iter < 100_000; iter++) {
            double res = idaSearch(g, start, 0, thresh, goal, onPath, pathEdges);
            if (res < 0)                   { found = true; break; }
            if (res == Double.MAX_VALUE)   break;
            thresh = res;
        }
        int useful = found ? pathEdges.size() + 1 : 0;
        double deltaD = 1.0;
        if (found && !pathEdges.isEmpty()) {
            double gCur = 0;
            for (int[] edge : pathEdges) {
                int u = edge[0], v = edge[1], ec = edge[2];
                List<int[]> nbrs = g.nbrs(u);
                double[] psiAll  = new double[nbrs.size()];
                for (int i = 0; i < nbrs.size(); i++)
                    psiAll[i] = gCur + nbrs.get(i)[2]
                              + g.h(g.id(nbrs.get(i)[0], nbrs.get(i)[1]));
                double gNext = gCur + ec;
                double[] a   = accuracyMeasures(psiAll, new double[]{gNext + g.h(v)});
                if (a[0] >= 0) deltaD *= a[0];
                gCur = gNext;
            }
        }
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = idaCount; m.useful = useful; m.found = found; m.par = null;
        m.deltaD = deltaD;
        m.deltaSMean = (idaStepsCount > 0) ? idaSumDS / idaStepsCount : 1.0;
        m.deltaCMean = 0.0;
        return m;
    }

    // ============================================================
    // ALGORITHM 9: ALPHA-BETA
    // ============================================================
    static long   abCount;
    static double abSumDS;
    static int    abSteps;

    static int alphaBeta(int[] tree, int node, int depth,
                          int alpha, int beta, boolean maxNode,
                          int b, int maxDepth, int firstLeaf) {
        abCount++;
        if (node >= firstLeaf || depth >= maxDepth)
            return tree[Math.min(node, tree.length - 1)];
        int[] children = childrenOf(node, b, tree.length);
        double[] psiAll = new double[children.length];
        for (int i = 0; i < children.length; i++)
            psiAll[i] = tree[Math.min(children[i], tree.length - 1)];
        List<Double> keptPsi = new ArrayList<>();

        if (maxNode) {
            Integer[] box = new Integer[children.length];
            for (int i = 0; i < children.length; i++) box[i] = children[i];
            Arrays.sort(box, (a2, b2) -> Integer.compare(
                    tree[Math.min(b2, tree.length-1)], tree[Math.min(a2, tree.length-1)]));
            for (int i = 0; i < children.length; i++) children[i] = box[i];
            int val = Integer.MIN_VALUE;
            for (int ch : children) {
                int cv = alphaBeta(tree, ch, depth+1, alpha, beta, false, b, maxDepth, firstLeaf);
                keptPsi.add((double) cv);
                val = Math.max(val, cv); alpha = Math.max(alpha, val);
                if (beta <= alpha) break;
            }
            double[] acc = accuracyMeasures(psiAll, keptPsi.stream().mapToDouble(x->x).toArray());
            if (acc[0] >= 0) { abSumDS += acc[0]; abSteps++; }
            return val;
        } else {
            Integer[] box = new Integer[children.length];
            for (int i = 0; i < children.length; i++) box[i] = children[i];
            Arrays.sort(box, (a2, b2) -> Integer.compare(
                    tree[Math.min(a2, tree.length-1)], tree[Math.min(b2, tree.length-1)]));
            for (int i = 0; i < children.length; i++) children[i] = box[i];
            int val = Integer.MAX_VALUE;            
            for (int ch : children) {
                int cv = alphaBeta(tree, ch, depth+1, alpha, beta, true, b, maxDepth, firstLeaf);
                keptPsi.add((double) cv);
                val = Math.min(val, cv);             
                beta = Math.min(beta, val);
                if (beta <= alpha) break;
            }
            double[] acc = accuracyMeasures(psiAll, keptPsi.stream().mapToDouble(x->x).toArray());
            if (acc[0] >= 0) { abSumDS += acc[0]; abSteps++; }
            return val;
        }
    }

    static int[] childrenOf(int node, int b, int treeLen) {
        int[] ch = new int[b];
        for (int i = 0; i < b; i++) {
            ch[i] = b * node + 1 + i;
            if (ch[i] >= treeLen) ch[i] = treeLen - 1;
        }
        return ch;
    }

    static EDBTMetrics abExpand(int b, int depth, long seed) {
        abCount = 0; abSumDS = 0.0; abSteps = 0;
        Random rng = new Random(seed);
        int safeDepth = Math.min(depth, AB_MAX_DEPTH);
        long size = 1, acc = 1;
        for (int i = 0; i < safeDepth; i++) {
            acc *= b;
            if (acc > 10_000_000L) { safeDepth = i; break; }
            size += acc;
        }
        int[] tree      = new int[(int) size];
        int   firstLeaf = (int) (size - acc);
        for (int i = firstLeaf; i < tree.length; i++) tree[i] = rng.nextInt(101);
        alphaBeta(tree, 0, 0, Integer.MIN_VALUE, Integer.MAX_VALUE,
                  true, b, safeDepth, firstLeaf);
        double deltaD = 1.0;
        int node = 0;
        for (int dep = 0; dep < safeDepth && node < firstLeaf; dep++) {
            int[]    ch     = childrenOf(node, b, tree.length);
            double[] psiAll = new double[ch.length];
            for (int i = 0; i < ch.length; i++) psiAll[i] = tree[Math.min(ch[i], tree.length-1)];
            double[] a = accuracyMeasures(psiAll, new double[]{psiAll[0]});
            if (a[0] >= 0) deltaD *= a[0];
            node = ch[0];
        }
        EDBTMetrics m = new EDBTMetrics();
        m.nodes = abCount; m.useful = depth + 1; m.found = true; m.par = null;
        m.deltaD = deltaD;
        m.deltaSMean = (abSteps > 0) ? abSumDS / abSteps : 1.0;
        m.deltaCMean = 0.0;
        return m;
    }

    // ============================================================
    // BETA ESTIMATION  
    // ============================================================
    static double[] fitBeta(List<Integer> depths, List<Double> meanNodes, int b) {
        int n = depths.size();
        if (n < 2) return new double[]{0.5, 1.0, 0.0};
        double logb = Math.log(b);
        double[] X = new double[n], Y = new double[n];
        for (int i = 0; i < n; i++) {
            X[i] = depths.get(i) * logb;
            Y[i] = Math.log(Math.max(meanNodes.get(i), 1.0));
        }
        double Xm = 0, Ym = 0;
        for (int i = 0; i < n; i++) { Xm += X[i]; Ym += Y[i]; }
        Xm /= n; Ym /= n;
        double sXX = 0, sXY = 0, sYY = 0;
        for (int i = 0; i < n; i++) {
            sXX += (X[i]-Xm)*(X[i]-Xm);
            sXY += (X[i]-Xm)*(Y[i]-Ym);
            sYY += (Y[i]-Ym)*(Y[i]-Ym);
        }
        double beta  = (sXX < 1e-12) ? 0.5 : sXY / sXX;
        beta         = clampBeta(beta);
        double a     = Ym - beta * Xm;
        double alpha = Math.exp(a);
        double r2    = (sYY < 1e-12) ? 1.0 : (sXY * sXY) / (sXX * sYY);
        return new double[]{beta, alpha, r2};
    }

    // ============================================================
    // EXPERIMENT RUNNER
    // ============================================================
    static List<AlgResult>           results   = new ArrayList<>();
    static Map<String,List<Double>>  algMeans  = new LinkedHashMap<>();
    static Map<String,List<Integer>> algDepths = new LinkedHashMap<>();
    static Map<String,List<Double>>  vizNodes  = new LinkedHashMap<>();

    static final String[] ALL_ALGS = {
        "BFS","Dijkstra","A*","WA*(1.5)","WA*(2.0)","GBFS",
        "BiDijkstra","RBFS","IDA*","AlphaBeta"
    };

    static void runExperiments() {
        System.out.println("=".repeat(80));
        System.out.println(" EDBT FRAMEWORK — EXPERIMENTAL VALIDATION");
        System.out.printf(" Grid=4d×4d | ρ=%.2f | ε=%.2f | Trials=%d%n%n",
                OBSTACLE_RATE, EPSILON, NUM_TRIALS);

        for (String a : ALL_ALGS) {
            algMeans.put(a, new ArrayList<>());
            algDepths.put(a, new ArrayList<>());
            vizNodes.put(a, new ArrayList<>());
        }

        for (int depth : DEPTHS) {
            int gs = 4 * depth;
            Map<String,Long>   sN  = new LinkedHashMap<>(), sN2 = new LinkedHashMap<>();
            Map<String,Long>   sU  = new LinkedHashMap<>(), sF  = new LinkedHashMap<>();
            Map<String,Double> sDD = new LinkedHashMap<>(), sDS = new LinkedHashMap<>(),
                               sDC = new LinkedHashMap<>();
            for (String a : ALL_ALGS) {
                sN.put(a,0L); sN2.put(a,0L); sU.put(a,0L); sF.put(a,0L);
                sDD.put(a,0.0); sDS.put(a,0.0); sDC.put(a,0.0);
            }

            for (int trial = 0; trial < NUM_TRIALS; trial++) {
                long seed = 100L * depth + trial;
                Grid g = new Grid(gs, seed);
                EDBTMetrics r;
                r = bfsExpand(g);               acc(sN,sN2,sU,sF,sDD,sDS,sDC,"BFS",r);
                r = dijkstraExpand(g);           acc(sN,sN2,sU,sF,sDD,sDS,sDC,"Dijkstra",r);
                r = astarExpand(g, 1.0);         acc(sN,sN2,sU,sF,sDD,sDS,sDC,"A*",r);
                r = astarExpand(g, 1.5);         acc(sN,sN2,sU,sF,sDD,sDS,sDC,"WA*(1.5)",r);
                r = astarExpand(g, 2.0);         acc(sN,sN2,sU,sF,sDD,sDS,sDC,"WA*(2.0)",r);
                r = gbfsExpand(g);               acc(sN,sN2,sU,sF,sDD,sDS,sDC,"GBFS",r);
                r = biDijkstraExpand(g);         acc(sN,sN2,sU,sF,sDD,sDS,sDC,"BiDijkstra",r);
                r = rbfsExpand(g);               acc(sN,sN2,sU,sF,sDD,sDS,sDC,"RBFS",r);
                r = idaExpand(g);                acc(sN,sN2,sU,sF,sDD,sDS,sDC,"IDA*",r);
                r = abExpand(BRANCHING,depth,seed);
                                                 acc(sN,sN2,sU,sF,sDD,sDS,sDC,"AlphaBeta",r);
            }

            System.out.printf("Depth=%2d (grid=%dx%d):%n", depth, gs, gs);
            for (String alg : ALL_ALGS) {
                double mn  = sN.get(alg)  / (double) NUM_TRIALS;
                double mu  = sU.get(alg)  / (double) NUM_TRIALS;
                double sr  = sF.get(alg)  / (double) NUM_TRIALS;
                double var = sN2.get(alg) / (double) NUM_TRIALS - mn * mn;
                double std = Math.sqrt(Math.max(var, 0));
                double dd  = sDD.get(alg) / (double) NUM_TRIALS;
                double ds  = sDS.get(alg) / (double) NUM_TRIALS;
                double dc  = sDC.get(alg) / (double) NUM_TRIALS;
                double q   = alg.startsWith("WA*(1.5)") ? 1.0/1.5
                           : alg.startsWith("WA*(2.0)") ? 1.0/2.0 : 1.0;
                double bei  = computeBEI(mn, mu, q, BRANCHING);
                double abei = computeABEI(mn, mu, q, dd, BRANCHING);

                List<Double>  mns = algMeans.get(alg);
                List<Integer> ds2 = algDepths.get(alg);
                mns.add(mn); ds2.add(depth);

                double[] fit  = fitBeta(ds2, mns, BRANCHING);
                double   beta = fit[0];
                double[] gam  = complexityProfile(beta, depth, BRANCHING, OBSTACLE_RATE);
                double gamP = gam[0], gamM = gam[1], gamMu = gam[2];
                double acp  = mn * dd;
                double cert = pcctCertificate(alg);
                double[] om = convergenceOperator(beta, dd, EPSILON, BRANCHING);

                AlgResult ar = new AlgResult(alg, depth, mn, std, mu, bei, sr,
                        dd, ds, dc, abei, om[1], cert, acp, gamM, gamP, gamMu);
                results.add(ar);
                vizNodes.get(alg).add(mn);

                System.out.printf(
                    "  %-14s N=%8.1f±%6.1f | Δᵈ=%.3f | δs=%.3f | " +
                    "BEI=%.4f | ABEI=%.4f | ACP=%.1f≥γ⁻=%.1f[%s] | q*=%.3f%n",
                    alg, mn, std, dd, ds, bei, abei, acp, gamM,
                    ar.acpSatisfied ? "✓" : "✗", om[1]);
            }
            System.out.println();
        }

        printBetaTable();
        printComplexityProfileTable();
        printACPTable();
        printPCCTTable();
        printABEIvsBEI();
        printConvergenceTable();
        printBEITable();
        printHeldOut();
        printRanking();
    }

    static void acc(Map<String,Long> sN, Map<String,Long> sN2,
                    Map<String,Long> sU, Map<String,Long> sF,
                    Map<String,Double> sDD, Map<String,Double> sDS,
                    Map<String,Double> sDC, String alg, EDBTMetrics m) {
        sN.merge(alg, m.nodes, Long::sum);
        sN2.merge(alg, m.nodes * m.nodes, Long::sum);
        sU.merge(alg, (long) m.useful, Long::sum);
        if (m.found) sF.merge(alg, 1L, Long::sum);
        sDD.merge(alg, m.deltaD, Double::sum);
        sDS.merge(alg, m.deltaSMean, Double::sum);
        sDC.merge(alg, m.deltaCMean, Double::sum);
    }

    // ============================================================
    // TABLES
    // ============================================================
    static void printBetaTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ESTIMATED β PARAMETERS");
        System.out.printf(" %-14s | %8s | %10s | %6s%n","Algorithm","β̂","α̂","R²");
        System.out.println(" "+"-".repeat(44));
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.isEmpty()) continue;
            double[] fit = fitBeta(d, mn, BRANCHING);
            System.out.printf(" %-14s | %8.4f | %10.4f | %6.4f%n", alg, fit[0], fit[1], fit[2]);
        }
        System.out.println();
    }

    static void printComplexityProfileTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: COMPLEXITY PROFILE BOUNDS");
        System.out.println("=".repeat(80));
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.isEmpty()) continue;
            double beta = fitBeta(d, mn, BRANCHING)[0];
            System.out.printf(" %s (β̂=%.4f):%n", alg, beta);
            System.out.printf("  %4s|%10s|%10s|%10s|%10s|%s%n","d","γ⁻","γᵘ","N̄","γ⁺","In[γ⁻,γ⁺]?");
            boolean allOk = true;
            for (AlgResult r : results) {
                if (!r.name.equals(alg)) continue;
                boolean ok = r.meanNodes >= r.gammaMinus - 1e-3
                          && r.meanNodes <= r.gammaPlus + r.gammaPlus * 0.05;
                if (!ok) allOk = false;
                System.out.printf("  %4d|%10.1f|%10.1f|%10.1f|%10.1f|%s%n",
                        r.depth, r.gammaMinus, r.gammaMu, r.meanNodes, r.gammaPlus, ok?"✓":"✗");
            }
            System.out.printf("  → All OK: %s%n%n", allOk ? "YES" : "NO");
        }
    }

    static void printACPTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ACP BOUND  N̄·Δᵈ ≥ γ⁻(β̂,d)");
        System.out.printf(" %-14s|%4s|%8s|%8s|%10s|%10s|%s%n",
                "Algorithm","d","N̄","Δᵈ","N̄·Δᵈ","γ⁻","OK?");
        System.out.println(" "+"-".repeat(68));
        for (AlgResult r : results)
            System.out.printf(" %-14s|%4d|%8.1f|%8.3f|%10.1f|%10.1f|%s%n",
                    r.name, r.depth, r.meanNodes, r.meanDeltaD,
                    r.acpProduct, r.gammaMinus, r.acpSatisfied?"✓":"✗");
        System.out.println();
    }

    static void printPCCTTable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: PCCT  δs_obs ≥ δs_cert");
        System.out.printf(" %-14s|%4s|%10s|%10s|%10s|%s%n",
                "Algorithm","d","δs_cert","δs_obs","Δᵈ","OK?");
        System.out.println(" "+"-".repeat(58));
        for (AlgResult r : results)
            System.out.printf(" %-14s|%4d|%10.3f|%10.3f|%10.3f|%s%n",
                    r.name, r.depth, r.certDeltaS, r.meanDeltaS,
                    r.meanDeltaD, r.pcctSatisfied?"✓":"✗");
        System.out.println();
    }

    static void printABEIvsBEI() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: ABEI vs BEI at d=10");
        System.out.printf(" %-14s|%8s|%8s|%8s|%8s|%s%n",
                "Algorithm","N̄","Δᵈ","BEI","ABEI","ABEI vs BEI");
        System.out.println(" "+"-".repeat(62));
        results.stream().filter(r -> r.depth == 10)
               .sorted(Comparator.comparingDouble((AlgResult r) -> r.meanABEI).reversed())
               .forEach(r -> System.out.printf(
                   " %-14s|%8.1f|%8.3f|%8.4f|%8.4f|%s%n",
                   r.name, r.meanNodes, r.meanDeltaD, r.meanBEI, r.meanABEI,
                   r.meanABEI > r.meanBEI ? "↑" : r.meanABEI < r.meanBEI ? "↓" : "="));
        System.out.println();
    }

    static void printConvergenceTable() {
        System.out.println("=".repeat(80));
        System.out.printf(" TABLE: CONVERGENCE OPERATOR Ω(ε=%.2f) at d=10%n", EPSILON);
        System.out.printf(" %-14s|%8s|%8s|%8s|%10s%n","Algorithm","β̂","Δᵈ","d*","q*");
        System.out.println(" "+"-".repeat(52));
        for (AlgResult r : results) {
            if (r.depth != 10) continue;
            List<Double> mn = algMeans.get(r.name); List<Integer> d = algDepths.get(r.name);
            double beta = fitBeta(d, mn, BRANCHING)[0];
            double[] om = convergenceOperator(beta, r.meanDeltaD, EPSILON, BRANCHING);
            System.out.printf(" %-14s|%8.4f|%8.3f|%8.1f|%10.4f%n",
                    r.name, beta, r.meanDeltaD, om[0], om[1]);
        }
        System.out.println();
    }

    static void printBEITable() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: BEI AND ABEI ACROSS DEPTHS");
        System.out.printf(" %-14s","Algorithm");
        for (int d : DEPTHS) System.out.printf("|d=%-2d BEI   ABEI ",d);
        System.out.println();
        System.out.println(" "+"-".repeat(14 + DEPTHS.length * 18));
        for (String alg : ALL_ALGS) {
            System.out.printf(" %-14s", alg);
            for (AlgResult r : results)
                if (r.name.equals(alg))
                    System.out.printf("|%6.4f %6.4f ", r.meanBEI, r.meanABEI);
            System.out.println();
        }
        System.out.println();
    }

    static void printHeldOut() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: HELD-OUT PREDICTION ERROR (train d=4,6,8; test d=10,12)");
        System.out.println("=".repeat(80));
        for (String alg : ALL_ALGS) {
            List<Double> mn = algMeans.get(alg); List<Integer> d = algDepths.get(alg);
            if (mn.size() < 5) continue;
            double[] fit     = fitBeta(d.subList(0,3), mn.subList(0,3), BRANCHING);
            double totalErr  = 0; int cnt = 0;
            for (int i = 3; i < d.size(); i++) {
                double pred = fit[1] * Math.pow(BRANCHING, fit[0] * d.get(i));
                double act  = mn.get(i);
                if (act > 0) { totalErr += Math.abs(pred - act) / act; cnt++; }
            }
            System.out.printf(" %-14s | MAPE = %6.2f%%%n",
                    alg, cnt > 0 ? totalErr / cnt * 100 : 0);
        }
        System.out.println();
    }

    static void printRanking() {
        System.out.println("=".repeat(80));
        System.out.println(" TABLE: RANKING AT d=12 (by ABEI)");
        System.out.println("=".repeat(80));
        results.stream().filter(r -> r.depth == 12)
               .sorted(Comparator.comparingDouble((AlgResult r) -> r.meanABEI).reversed())
               .forEach(new java.util.function.Consumer<AlgResult>() {
                   int rank = 1;
                   public void accept(AlgResult r) {
                       System.out.printf(
                           " #%d %-14s N=%9.1f | BEI=%.4f | ABEI=%.4f | Δᵈ=%.3f | q*=%.3f%n",
                           rank++, r.name, r.meanNodes, r.meanBEI,
                           r.meanABEI, r.meanDeltaD, r.meanQstar);
                   }
               });
        System.out.println();
    }

    // ============================================================
    // VISUALIZATION
    // ============================================================
    static final Color[] ALG_COLORS = {
        new Color(220,50,50),  new Color(160,32,240), new Color(30,144,255),
        new Color(0,180,80),   new Color(50,200,100), new Color(255,165,0),
        new Color(0,206,209),  new Color(255,20,147), new Color(210,180,140),
        new Color(139,90,43)
    };

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        int W = getWidth(), H = getHeight();
        g2.setColor(new Color(250,250,250)); g2.fillRect(0,0,W,H);
        int lW = (int)(W * 0.62);
        drawLineChart(g2, 0, 0, lW, H);
        drawBEIBars(g2, lW, 0, W - lW, H);
    }

    void drawLineChart(Graphics2D g2, int ox, int oy, int W, int H) {
        int mL=70,mR=20,mT=50,mB=55,pW=W-mL-mR,pH=H-mT-mB;
        g2.setColor(Color.WHITE); g2.fillRect(ox+mL,oy+mT,pW,pH);
        g2.setColor(Color.BLACK);
        g2.setFont(new Font("SansSerif",Font.BOLD,12));
        g2.drawString("Node Expansions vs Depth (log scale)",ox+mL,oy+30);
        double logMax=1,logMin=0;
        for (List<Double> v : vizNodes.values())
            for (double x : v) if (x>1){
                logMax=Math.max(logMax,Math.log10(x));
                logMin=Math.min(logMin,Math.log10(x));
            }
        logMax=Math.ceil(logMax)+0.5; logMin=Math.max(0,Math.floor(logMin)-0.5);
        double lr=logMax-logMin;
        g2.setFont(new Font("SansSerif",Font.PLAIN,9));
        for (int pw=(int)logMin;pw<=(int)logMax;pw++) {
            int yy=oy+mT+pH-(int)(((pw-logMin)/lr)*pH);
            if (yy<oy+mT||yy>oy+mT+pH) continue;
            g2.setColor(new Color(210,210,210));
            g2.setStroke(new BasicStroke(1,BasicStroke.CAP_BUTT,BasicStroke.JOIN_BEVEL,
                    0,new float[]{4},0));
            g2.drawLine(ox+mL,yy,ox+mL+pW,yy);
            g2.setColor(Color.GRAY); g2.setStroke(new BasicStroke(1));
            g2.drawString("10^"+pw,ox+mL-42,yy+4);
        }
        g2.setColor(Color.BLACK); g2.setStroke(new BasicStroke(2));
        g2.drawLine(ox+mL,oy+mT,ox+mL,oy+mT+pH);
        g2.drawLine(ox+mL,oy+mT+pH,ox+mL+pW,oy+mT+pH);
        int nd=DEPTHS.length;
        g2.setFont(new Font("SansSerif",Font.PLAIN,10));
        for (int i=0;i<nd;i++) {
            int x=ox+mL+(i*pW)/(nd-1);
            g2.setStroke(new BasicStroke(1));
            g2.drawLine(x,oy+mT+pH,x,oy+mT+pH+4);
            g2.drawString("d="+DEPTHS[i],x-10,oy+mT+pH+16);
        }
        g2.setFont(new Font("SansSerif",Font.BOLD,10));
        g2.drawString("Search Depth",ox+mL+pW/2-35,oy+mT+pH+35);
        int ai=0;
        for (String alg : ALL_ALGS) {
            List<Double> vals=vizNodes.get(alg);
            if (vals==null||vals.isEmpty()){ai++;continue;}
            g2.setColor(ALG_COLORS[ai]); g2.setStroke(new BasicStroke(2));
            int px=-1,py=-1;
            for (int i=0;i<vals.size();i++) {
                double v=vals.get(i); if (v<1) continue;
                int x=ox+mL+(i*pW)/(nd-1);
                int y=oy+mT+pH-(int)(((Math.log10(v)-logMin)/lr)*pH);
                y=Math.max(oy+mT,Math.min(oy+mT+pH,y));
                g2.fillOval(x-4,y-4,8,8);
                if (px>=0) g2.drawLine(px,py,x,y);
                px=x;py=y;
            }
            ai++;
        }
        int lx=ox+mL+10,ly=oy+mT+10;
        g2.setColor(new Color(255,255,255,200));
        g2.fillRoundRect(lx-4,ly-4,135,ALL_ALGS.length*16+8,6,6);
        for (int i=0;i<ALL_ALGS.length;i++) {
            g2.setColor(ALG_COLORS[i]); g2.setStroke(new BasicStroke(2));
            g2.drawLine(lx,ly+i*16+6,lx+18,ly+i*16+6);
            g2.fillOval(lx+5,ly+i*16+2,8,8);
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("SansSerif",Font.PLAIN,9));
            g2.drawString(ALL_ALGS[i],lx+22,ly+i*16+10);
        }
    }

    void drawBEIBars(Graphics2D g2, int ox, int oy, int W, int H) {
        int mL=85,mR=15,mT=50,mB=55,pW=W-mL-mR,pH=H-mT-mB;
        g2.setColor(Color.WHITE); g2.fillRect(ox,oy,W,H);
        g2.setColor(Color.BLACK);
        g2.setFont(new Font("SansSerif",Font.BOLD,12));
        g2.drawString("ABEI per Algorithm & Depth",ox+mL,oy+30);
        double maxB=results.stream().mapToDouble(r->r.meanABEI).max().orElse(0.1);
        if (maxB<1e-9) maxB=0.1;
        int nD=DEPTHS.length,nA=ALL_ALGS.length;
        int gW=pW/nD,bW=Math.max(2,(gW-6)/nA);
        g2.setColor(Color.BLACK); g2.setStroke(new BasicStroke(2));
        g2.drawLine(ox+mL,oy+mT,ox+mL,oy+mT+pH);
        g2.drawLine(ox+mL,oy+mT+pH,ox+mL+pW,oy+mT+pH);
        g2.setFont(new Font("SansSerif",Font.PLAIN,8));
        for (int t=0;t<=5;t++) {
            double val=maxB*t/5.0; int yy=oy+mT+pH-(int)(pH*t/5.0);
            g2.setColor(new Color(200,200,200));
            g2.setStroke(new BasicStroke(1,BasicStroke.CAP_BUTT,BasicStroke.JOIN_BEVEL,
                    0,new float[]{3},0));
            g2.drawLine(ox+mL,yy,ox+mL+pW,yy);
            g2.setColor(Color.GRAY); g2.setStroke(new BasicStroke(1));
            g2.drawString(String.format("%.3f",val),ox+mL-50,yy+3);
        }
        for (int di=0;di<nD;di++) {
            int gx=ox+mL+di*gW+3, d=DEPTHS[di];
            for (int ai=0;ai<nA;ai++) {
                String alg=ALL_ALGS[ai]; double abei=0;
                for (AlgResult r:results) if (r.name.equals(alg)&&r.depth==d) abei=r.meanABEI;
                int bh=(int)(pH*abei/maxB),bx=gx+ai*bW,by=oy+mT+pH-bh;
                g2.setColor(ALG_COLORS[ai]); g2.fillRect(bx,by,Math.max(bW-1,1),bh);
                g2.setColor(ALG_COLORS[ai].darker()); g2.setStroke(new BasicStroke(1));
                g2.drawRect(bx,by,Math.max(bW-1,1),bh);
            }
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("SansSerif",Font.PLAIN,9));
            g2.drawString("d="+d,gx+gW/2-10,oy+mT+pH+14);
        }
        g2.setFont(new Font("SansSerif",Font.BOLD,10));
        g2.drawString("Search Depth",ox+mL+pW/2-35,oy+mT+pH+35);
    }

    // ============================================================
    // MAIN
    // ============================================================
    public static void main(String[] args) {
        System.out.println("EDBT 8-tuple: T⁺=(V,E,Φ,Ψ,Π,γ⃗,Δ,Ω)");
        runExperiments();
        SwingUtilities.invokeLater(() -> {
            JFrame f = new JFrame("EDBT Framework Results");
            f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            f.setSize(1400,700);
            UnifiedTreeExperimentL p = new UnifiedTreeExperimentL();
            p.setBackground(Color.WHITE);
            f.add(p);
            f.setLocationRelativeTo(null);
            f.setVisible(true);
        });
    }
}