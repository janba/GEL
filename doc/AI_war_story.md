# Some anecdotes about the use of generative AI in Programming

## What AI?

The AI used was Grok 4.7 (and earlier 4.5 and 4.6) used via the command line interface. This is supposedly less powerful than the frontier models from Anthropic and OpenAI but not that far behind.

## How it works:

I use Grok CLI which runs in the terminal. Being curious about exactly where things happen, I asked Grok. The split is:

• The model runs on xAI's servers. Your message is sent there, and the reply, including any request to run a command or edit a file, comes back from there.
• The Grok app on this computer carries out those requests. Shell commands, file reads, and file edits happen locally in /Users/janba/GEL and under your home directory. The results are sent back so the model can continue.
• A few tools are remote even when the decision to call them is made from here: web search, fetching a URL, X search, and image generation. The LRsR runs did not use those.

## Grok vs Github Copilot. Round 1, proofreading.

I have updated the book with the changes from Springer's proofreading process. I had Github Copilot go over the word files and the tex files and find the proofreading changes and then prepare a script to apply the corrections. I also asked it to find additional issues and correct those. It found relatively few things.

This was surprising since Jens and I had done a test ono Chapter 1 this Friday using Grok which seemed to find quite a lot of changes just in Chapter 1.

So I had Grok redo the comparison between word and tex, and it found many more things to change. It also found more changes that were not in the word files.

I then asked both AIs why they had a different result. They don't quite agree.

Copilot: The practical difference is: my work was “evidence-based, minimal-risk corrections,” while the Grok run looks more like “comprehensive editorial proofreading.”

Grok: Yes. Copilot did not really do the three tasks. It did a shallow typo grep on a few files and stopped.

I think this is so funny. Copilot is so smooth and a master of spin. Grok is more passive aggressive but also more diligent. 

## Grok vs Github Copilot. Round 2, Curvature estimation

I was changing the curvature tensor estimation from only finding a paraboloid with no linear terms to first fitting a plane to the local neighborhood and then a paraboloid. This is a bit of a hack, but it should work about as well as fitting a paraboloid with linear terms. In practice it did not, so I asked both Grok and Copilot about that.

Copilot: The fix is to estimate the linear and quadratic terms together, with 5 unknowns (x, y, ½x², xy, ½y²) in a single least-squares fit.
... Copilot said much more but that was the recommendation. It all take 19 seconds. The model was Claude Sonnet 5.5 and no files other than the one with the problematic functions were consulted.

Grok thought for 16 minutes and 27 seconds and came back with largely the same recommendation. However, it also found a bug in an unrelated downstream piece of code. The way I compute eigenvalues (using power iterations) does not work reliably when two eigenvalues are numerically the same but of opposite sign. It came up with a solution for that problem as well (shifting such that the eigenvalues are different) but also told me that the eigenvalues should be computed directly since it is a 2x2 matrix. 

So, Grok did much more work again. Arguably, though, Copilot was more efficient.

## Build System

AI has been used to make several changes to the way GEL is built:
- Make a wheel per platform rather than a single wheel which loads the right binary.
- automatic deployment to PyPI. When the wheels are built using Github actions they are also uploaded to PyPI if the commit is tagged with a version number.
- More general streamlining - turning GEL into a CMake package. Make demos build in the same way that an ordinary project using GEL should be built.
- Produce a yaml file which tests installation on a range of platforms and runs a simple test of GEL features.

Interestingly the install test passed despite the AI having initially produced a wheel that would only really install on computers with almost identical setup to mine. That was because the installation just installed a previous version.

#  Developing a reconstruction method

Some years ago we developed the RsR method for combinatorial triangle mesh reconstruction from point clouds. While it works really well, it has a couple of drawbacks. First of all, it is not very fast. Second, it is reliant on consistently oriented normals. We tried to address the first drawback with a hierarchical method. As of writing this method is unpublished, and while it is somewhat faster it is not dramatically faster and there is some cost to quality. Regarding the second issue, I first tried to fix the normals as discussed below.

The following is a brief account of how I used generative AI to improve the method and was gently led towards implementing a slightly different method. I will try to describe how AI was used and how it influenced the work process. 

## First Analysis

The first step was the recognition that RsR relied on point normals and our method for computing point normals was not good enough. A point normal is a vector that is estimated to be perpendicular to the local tangent plane. Typically, the tangent plane is estimated via PCA. This process is effective, but it relies on a neighborhood size. Moreover, it finds the tangent plane normal but not the orientation. This has to be done in a subsequent step. Our method for computing normals was slightly naive and worked poorly when the point cloud was sparse. Moreover, the orientation computation was brittle.

## Improving Normal Computation
Thus, the first order of business was to improve normal computation. Various strategies for computing normals were discussed with Grok via the web chat interface. Having considered more advanced methods, we finally settled on a traditional method, but with a twist. PCA is run on several neighborhood sizes and the entropy is computed for each. Entropy is minimized when the smallest eigenvector of the covariance matrix is much smaller than the two others. I realized that this worked poorly in some cases: when the points are on a line rather than distributed over a disk-like patch. There was already a fix in place but it only worked in extreme cases. Grok proposed a method for making the entropy higher for the close to degenerate cases than the ones where the two large eigenvectors are similar in size. With this in place, we had a fairly good normal estimaton which adapted to the local feature size. Sadly this did not fix the problem that consistent normal orientation computation was brittle. Now, my PhD student Ruiqi had (I believe) already implemented a distance weighting that assigns larger distances to edges whose vertex normals are not very perpendicular to the edges. This is used when a minimum spanning tree is computed for the points in a point cloud. It assigns longer distances to pairs of points that are unlikely to have similar normals which means that such edges are not in the MST. That is useful when we propagate the normal orientation along the MST.

## Researching Delaunay and Gabriel Graphs

Now, all problems with normal computation were sort of solved, but the original RsR is still a bit brittle in a few cases where the input is a set of points from a triangle mesh, and hRsR felt decidedly unsatisfactory.

This made me want to refactor the method. I started by asking Grok to implement a very efficient method for 3D Delaunay graph computation that Bruno Levy has (probably) invented and included in Geogram. The method is relatively simple to implement and I described it to Grok who implemented it with few issues. It also implemented a similar method for finding Gabriel edges (edges that are diameters in their Delaunay circumsphere) and it improved the method for minimum spanning tree computation. Finally, Grok made a simple viewer based on existing PyGEL components for inspecting the resulting graphs and trees.

This was very helpful for human understanding of the problem. I could now look at the Delaunay graph, Gabriel graph and MST and get a sense of what algorithms might work.

## Improving the kD Tree

While thinking about what to do, I had Grok improve the kD tree. The key thing it did was to change the tree from a binary heap to a more traditional tree structure with at most 16 points in the leaf nodes. This improved efficiency of both building and lookup by around 20% at some cost in terms of memory.

When I tested this it seemed to be slower. Grok found out that I was running a debug build. This was surprising since the CMake files always make Release builds unless otherwise instructed. Grok figured out that an old CMake run had produced Debug make files and since a new CMake run had not been executed we still had the Debug setting. This is now fixed by having the install script insist on deciding whether to build debug or release. I think wrangling build files is one of the areas where the AI help is even more welcome than in other places.

In some cases, the code generated was quite sophisticated. In the code below, a lambda function is used to fill data into a container in order to allow for multiple interfaces.

```
unsigned in_sphere(const KeyType& p,
                ScalarType dist,
                std::vector<KeyT>& keys,
                std::vector<ValT>& vals) const
{
    assert(is_built);
    if(nodes.size()>1)
    {
        keys.clear();
        vals.clear();
        ScalarType max_sq_dist = CGLA::sqr(dist);
        auto sink = [&](const KDNode& node) {
            keys.push_back(node.key);
            vals.push_back(node.val);
        };
        in_sphere_priv(1, p, max_sq_dist, sink);
        return static_cast<unsigned>(keys.size());
    }
    return 0;
}
```

## LRsR: Attempt 1

Now, I was ready to make a new program for surface reconstruction - or so I thought. The idea was to wake the minimum spanning tree and create a triangle fan at each vertex. The triangles would be taken from the Delaunay triangulation which Grok had already coded up. of course, the triangles should be consistent when two one rings overlapped.

This method worked really poorly. In retrospect, I had not thought it completely through and there were too many places where it was not completely clear what to do. It is worth remarking that in such cases, AI in general and Grok in particular will cheerfully implement the closest approximation of what you ask regardless of whether it makes the slightest sense or not. After a few attempts at making this code work reasonably, I simply gave up. It is conceivable the method could have been made to work, but I had a better idea at that point.

## LRsR: Attempt 2 (less rotation based)

The idea was to start from the edges of the minimum spanning tree. Those edges should definitely be in the output mesh. If we could find good triangular faces incident on the edges of the tree that would be a start. This turned out to work really well. Now there is a lot of conditions that you can add to ensure you get a sane mesh. You want only two faces incident on each edge and you want the triangles to form a single fan around each vertex. Finally, it is easy to assign a quality value to triangle faces. Combining these things I arrived at a good selection of faces - and importantly a mostly consistent selection of faces. However, there were a lot of gaps since this method would only find faces incident on the minimum spanning tree edges.

Then I had the idea to ask Grok to add the newly generated edges, i.e. the triangle edges not in the MST, to the list of edges from which to find triangles. With this in place, the algorithm would find a quite good set of triangles. However, there were some defects. In its state at that time, the algorithm found triangles in no particular order. For a given edge a quality measure was used to find the best triangles but this quality measure was not used for global prioritization. The next step was to ask Grok to make a prority queue of triangles such that they would be added in descending order of quality. This was a straight forward change. With this in place, I got a pretty good triangulation of the point set I used for testing. This point set was quite challenging and contained very thin structures, but it now worked except that the triangles did not have a consistent winding order. However, since the triangles were all edge connected, it was easy for Grok to make code for consistently orienting the triangles.  With this in place a good mesh was produced. It still contains a few holes, but these are easy to close.

## GGSR

The Gabriel graph is very fast to compute. For each edge we only need to check whether it is the diameter of a sphere that contains no points. This is much faster than finding Delaunay edges, so if we could perform reconstruction just with Gabriel edges we are golden. 

A new prototype was made. In the first step, the gabriel graph was found and we tried to find all cycles in the mesh. At first, I described a homegrown algorithm for finding cycles based on running Dijkstra at each vertex. Grok implemented it, and it kind of worked but not very efficiently.

Using a chatbot for research I arrived at Johnson 1975 as the most obvious candidate algorithm. This was duly implemented and worked better but still the cycle finding made the whole thing much slower than the Delaunay graph based method from the first prototype.

Now, Johnson 1975 is a generic algorithm that does not exploit the geometry of the problem.

I described another approach to Grok which was based on first finding corners (a vertex and two relatively close incident edges) and then cycles by tracing a loop of corners where the last edge of one corner is the first edge of the next corner.

This was much more efficient, but ultimately some of the desired faces from the Gabriel Graph remain very hard to find.

## Main takeaways:
- AI is a godsend when it comes to build files. Don't touch them, just let AI deal with that. However, be aware that what you ask for might not be exactly what you want.
- Generative AI makes the entire slice - right from the kernel changes to the UI. That is super helpful.
- I don't really look at the code except maybe to make one-line changes. I think we have to trust the AI but validate the outcome. 
- The AI does not come up with the strategy for solving the problem. It is useful for researching various methods to achieve some goal, but I find it less useful for ultimately choosing the right strategy. This could be very problem-dependent and it could depend on the human and on the AI, but I think that as long as humans formulate the problem, we might have the sharpest intuition ultimately for the solution.
- If I don't understand what I am asking, the AI can't code it. My first solution did not work and Grok could not fix my idea. I have seen this before when I tried to use AI to generate ML-based texture synthesis (mostly Copilot though). The programs would work, but they would not synthesize textures.
- Since you do not author the code, throwing it away becomes much less painful. This is perhaps the benefit as it leaves you free to explore. Ultimately, the Delaunay Graph is a better starting point for mesh generation it seems now and the ability to play with the algorithms really helped bring clarity.
