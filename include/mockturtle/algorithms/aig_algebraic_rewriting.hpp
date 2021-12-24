/*!
  \file aig_algebraic_rewriting.hpp
  \brief AIG algebraric rewriting

  EPFL CS-472 2021 Final Project Option 1
*/

#pragma once

#include "../networks/aig.hpp"
#include "../views/depth_view.hpp"
#include "../views/topo_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class aig_algebraic_rewriting_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  aig_algebraic_rewriting_impl( Ntk& ntk )
    : ntk( ntk )
  {
    static_assert( has_level_v<Ntk>, "Ntk does not implement depth interface." );
  }

  void run()
  {
    bool cont{true}; /* continue trying */
    while ( cont )
    {
      cont = false; /* break the loop if no updates can be made */
      ntk.foreach_gate( [&]( node n ){
        if ( try_algebraic_rules( n ) )
        {
          ntk.update_levels();
          cont = true;
        }
      });
    }
  }

private:
  /* Try various algebraic rules on node n. Return true if the network is updated. */
  bool try_algebraic_rules( node n )
  {
    if ( try_associativity( n ) )
      return true;
    if ( try_distributivity( n ) )
      return true;
    /* TODO: add more rules here... */
    if ( try_three_level_distributivity( n ) )
      return true; 
    return false;
  }

  /* Try the associativity rule on node n. Return true if the network is updated. */
  bool try_associativity( node n )
  {
    /* TODO */
	if(ntk.fanin_size(n) && ntk.is_on_critical_path(n))
	{
		std::vector<signal> n_fanin_node_signals;
		node a, b;
		//First we have to check if a fanin of n is on the critical path.
		ntk.foreach_fanin(n, [&](auto const f)
		{
			n_fanin_node_signals.push_back(f);
		});
		a = ntk.get_node(n_fanin_node_signals[0]);
		b = ntk.get_node(n_fanin_node_signals[1]);

    		if(!((ntk.is_pi(a)) ^ (ntk.is_pi(b))))  
		{
      			return false;
    		}

		//if the node on the critical path is complemented, swapping cannot be done because we have a(cd)' and we cannot decompose (cd)'.
		if(ntk.is_on_critical_path(a) && !(ntk.is_complemented(n_fanin_node_signals[0])))
		{
			if(ntk.fanin_size(a))
			{
				//checking if any fanins of a are on the critical path.
				std::vector<signal> a_fanin_node_signals;
				node c, d;
				ntk.foreach_fanin(a, [&](auto const f)
				{
					a_fanin_node_signals.push_back(f);
				});
				c = ntk.get_node(a_fanin_node_signals[0]);
				d = ntk.get_node(a_fanin_node_signals[1]);

				if(ntk.is_on_critical_path(c) && !(ntk.is_on_critical_path(d)) && ((ntk.level(c) - ntk.level(b)) >= 1) && !(ntk.is_pi(c)))
		                {
					//c is on critical path.
					//and(and(b, d), c)
					auto aig1 = ntk.create_and(n_fanin_node_signals[1], a_fanin_node_signals[1]);
					auto aig2 = ntk.create_and(a_fanin_node_signals[0], aig1);
					ntk.substitute_node(n, aig2);
					return true;
				}
				else if(ntk.is_on_critical_path(d) && !(ntk.is_on_critical_path(c)) && ((ntk.level(d) - ntk.level(b)) >= 1) && !(ntk.is_pi(d)))
		                {
					//d is on critical path
					//and(and(b, c), d)
					auto aig1 = ntk.create_and(n_fanin_node_signals[1], a_fanin_node_signals[0]);
                                        auto aig2 = ntk.create_and(a_fanin_node_signals[1], aig1);
                                        ntk.substitute_node(n, aig2);
                                        return true;
				}
				//else do nothing.

			}
		}
		else if(ntk.is_on_critical_path(b) && !(ntk.is_complemented(n_fanin_node_signals[1])))
		{
			if(ntk.fanin_size(b))
			{
				//checking if any fanins of a are on the critical path.
				std::vector<signal> b_fanin_node_signals;
                                node c, d;
                                ntk.foreach_fanin(b, [&](auto const f)
                                {
                                        b_fanin_node_signals.push_back(f);
                                });
                                c = ntk.get_node(b_fanin_node_signals[0]);
                                d = ntk.get_node(b_fanin_node_signals[1]);

                                if(ntk.is_on_critical_path(c) && !(ntk.is_on_critical_path(d)) && ((ntk.level(c) - ntk.level(a)) >= 1) && !(ntk.is_pi(c)))
                                {
                                        //c is on critical path.
                                        //and(and(a, d), c)
                                        auto aig1 = ntk.create_and(n_fanin_node_signals[0], b_fanin_node_signals[1]);
                                        auto aig2 = ntk.create_and(b_fanin_node_signals[0], aig1);
                                        ntk.substitute_node(n, aig2);
                                        return true;
                                }
                                else if(ntk.is_on_critical_path(d) && !(ntk.is_on_critical_path(c)) && ((ntk.level(d) - ntk.level(a)) >= 1) && !(ntk.is_pi(d)))
                                {
                                        //d is on critical path
                                        //and(and(a, c), d)
                                        auto aig1 = ntk.create_and(n_fanin_node_signals[0], b_fanin_node_signals[0]);
                                        auto aig2 = ntk.create_and(b_fanin_node_signals[1], aig1);
                                        ntk.substitute_node(n, aig2);
                                        return true;
                                }
                                //else do nothing.
			}
		}
		//else do nothing.		

	}  //check that node has fanins and is not a constant or ci.

	return false;
  
  }

  /* Try the distributivity rule on node n. Return true if the network is updated. */
  bool try_distributivity( node n )
  {
    /* TODO */
        if(ntk.fanin_size(n) && ntk.is_on_critical_path(n))
        {

                std::vector<signal> n_fanin_node_signals;
                node a, b;
                //First we have to check if a fanin of n is on the critical path.
                ntk.foreach_fanin(n, [&](auto const f)
                {
                        n_fanin_node_signals.push_back(f);
                });
                a = ntk.get_node(n_fanin_node_signals[0]);
                b = ntk.get_node(n_fanin_node_signals[1]);

                if((ntk.is_pi(a)) || (ntk.is_pi(b)))
                {
                        return false;
                }

		//Step 1: To detect an OR (using complements) & Step 2: Find an intersection in fanins of fanins of n.
		if(ntk.is_complemented(n_fanin_node_signals[0]) && ntk.is_complemented(n_fanin_node_signals[1]) && ntk.fanin_size(a) && ntk.fanin_size(b))
		{
			std::vector<signal> a_fanin_node_signals;
			node c, d;
			ntk.foreach_fanin(a, [&](auto const f)
			{
				a_fanin_node_signals.push_back(f);
			});
			c = ntk.get_node(a_fanin_node_signals[0]);
			d = ntk.get_node(a_fanin_node_signals[1]);

			std::vector<signal> b_fanin_node_signals;
                        node x, y;
                        ntk.foreach_fanin(b, [&](auto const f)
                        {
                                b_fanin_node_signals.push_back(f);
                        });
                        x = ntk.get_node(b_fanin_node_signals[0]);
                        y = ntk.get_node(b_fanin_node_signals[1]);


			if(c == x && ntk.is_complemented(a_fanin_node_signals[0]) == ntk.is_complemented(b_fanin_node_signals[0]) && ntk.is_on_critical_path(c))	//c & x
			{
				auto aig1 = ntk.create_and(!a_fanin_node_signals[1], !b_fanin_node_signals[1]);
				auto aig2 = ntk.create_and(a_fanin_node_signals[0], !aig1);
				ntk.substitute_node(n, !aig2);
				return true;

			}
			else if(d == x && ntk.is_complemented(a_fanin_node_signals[1]) == ntk.is_complemented(b_fanin_node_signals[0]) && ntk.is_on_critical_path(d)) 	//d & x
			{
				auto aig1 = ntk.create_and(!a_fanin_node_signals[0], !b_fanin_node_signals[1]);
                                auto aig2 = ntk.create_and(a_fanin_node_signals[1], !aig1);
                                ntk.substitute_node(n, !aig2);
                                return true;
			}
			else if(c == y && ntk.is_complemented(a_fanin_node_signals[0]) == ntk.is_complemented(b_fanin_node_signals[1]) && ntk.is_on_critical_path(c))	//c & y
			{
				auto aig1 = ntk.create_and(!a_fanin_node_signals[1], !b_fanin_node_signals[0]);
                                auto aig2 = ntk.create_and(a_fanin_node_signals[0], !aig1);
                                ntk.substitute_node(n, !aig2);
                                return true;
			}
			else if(d == y && ntk.is_complemented(a_fanin_node_signals[1]) == ntk.is_complemented(b_fanin_node_signals[1]) && ntk.is_on_critical_path(d))	//d & y
			{
				auto aig1 = ntk.create_and(!a_fanin_node_signals[0], !b_fanin_node_signals[0]);
                                auto aig2 = ntk.create_and(a_fanin_node_signals[1], !aig1);
                                ntk.substitute_node(n, !aig2);
                                return true;
			}
			//else no intersection
		}
	}
	return false;
  }

  bool try_three_level_distributivity( node n )
  {
	if(ntk.fanin_size(n) && ntk.is_on_critical_path(n))
        {

                std::vector<signal> n_fanin_node_signals;
                node a, b;
                //First we have to check if a fanin of n is on the critical path.
                ntk.foreach_fanin(n, [&](auto const f)
                {
                        n_fanin_node_signals.push_back(f);
                });
                a = ntk.get_node(n_fanin_node_signals[0]);
                b = ntk.get_node(n_fanin_node_signals[1]);

                if((ntk.is_pi(a)) && (ntk.is_pi(b)))
                {
                        return false;
                }

		/* This rule is not in the pdf, but also simple:
     		((g x2) + x3 ) x4 = (g x2 x4) + (x3 x4) = (g (x2 x4)) + (x3 x4) */
		
		//Step 1: Identify critical path fanin viz. also complemented 
		//Step 2: move further into the critical path depth, and confirm that the level of x4 increased by 2 is still less than initial level of g.
		//Step 3: Apply the optimization

		if(ntk.is_on_critical_path(a) && ntk.is_complemented(n_fanin_node_signals[0]) && !ntk.is_on_critical_path(b))
		{
			//a is on critical path, b = x4. 
			std::vector<signal> a_fanin_node_signals;
                        node c, d;
                        ntk.foreach_fanin(a, [&](auto const f)
                        {
                                a_fanin_node_signals.push_back(f);
                        });
                        c = ntk.get_node(a_fanin_node_signals[0]);
                        d = ntk.get_node(a_fanin_node_signals[1]);

			if(ntk.is_complemented(a_fanin_node_signals[0]) && ntk.is_complemented(a_fanin_node_signals[1]))
			{
				if(ntk.is_on_critical_path(c) && !ntk.is_on_critical_path(d))
				{
					//d = x3
					std::vector<signal> c_fanin_node_signals;
					node x, y;
					ntk.foreach_fanin(c, [&](auto const f)
					{
						c_fanin_node_signals.push_back(f);
					});
					x = ntk.get_node(c_fanin_node_signals[0]);
					y = ntk.get_node(c_fanin_node_signals[1]);

					if(ntk.is_on_critical_path(x) && !ntk.is_on_critical_path(y))
					{
						//y = x2, x = g
						//if(ntk.level(x) > (ntk.level(b) + 2))
						{
							auto aig1 = ntk.create_and(c_fanin_node_signals[1], n_fanin_node_signals[1]); //x2 x4
							auto aig2 = ntk.create_and(c_fanin_node_signals[0], aig1); // g (x2 x4)
							auto aig3 = ntk.create_and(!a_fanin_node_signals[1], n_fanin_node_signals[1]); //x3 x4
							auto aig4 = ntk.create_and(!aig2, !aig3); 

							ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
							return true;
						}
						//else no benefit in terms of depth	
					}
					else if(ntk.is_on_critical_path(y) && !ntk.is_on_critical_path(x))
					{
						//x = x2, y = g
						//if(ntk.level(y) > (ntk.level(b) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(c_fanin_node_signals[0], n_fanin_node_signals[1]); //x2 x4
                                                        auto aig2 = ntk.create_and(c_fanin_node_signals[1], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!a_fanin_node_signals[1], n_fanin_node_signals[1]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
					}
				}
				else if(ntk.is_on_critical_path(d) && !ntk.is_on_critical_path(c))
				{
					//c = x3
					std::vector<signal> d_fanin_node_signals;
                                        node x, y;
                                        ntk.foreach_fanin(d, [&](auto const f)
                                        {
                                                d_fanin_node_signals.push_back(f);
                                        });
                                        x = ntk.get_node(d_fanin_node_signals[0]);
                                        y = ntk.get_node(d_fanin_node_signals[1]);

                                        if(ntk.is_on_critical_path(x) && !ntk.is_on_critical_path(y))
                                        {
                                                //y = x2, x = g
                                                //if(ntk.level(x) > (ntk.level(b) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(d_fanin_node_signals[1], n_fanin_node_signals[1]); //x2 x4
                                                        auto aig2 = ntk.create_and(d_fanin_node_signals[0], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!a_fanin_node_signals[0], n_fanin_node_signals[1]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
                                        }
                                        else if(ntk.is_on_critical_path(y) && !ntk.is_on_critical_path(x))
                                        {
                                                //x = x2, y = g
                                                //if(ntk.level(y) > (ntk.level(b) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(d_fanin_node_signals[0], n_fanin_node_signals[1]); //x2 x4
                                                        auto aig2 = ntk.create_and(d_fanin_node_signals[1], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!a_fanin_node_signals[0], n_fanin_node_signals[1]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
                                        }
				}
			}
			//else no distribution rule
		}
		else if(ntk.is_on_critical_path(b) && ntk.is_complemented(n_fanin_node_signals[1]) && !ntk.is_on_critical_path(a))
                {
                        //b is on critical path, a = x4.
			std::vector<signal> b_fanin_node_signals;
                        node c, d;
                        ntk.foreach_fanin(b, [&](auto const f)
                        {
                                b_fanin_node_signals.push_back(f);
                        });
                        c = ntk.get_node(b_fanin_node_signals[0]);
                        d = ntk.get_node(b_fanin_node_signals[1]);

                        if(ntk.is_complemented(b_fanin_node_signals[0]) && ntk.is_complemented(b_fanin_node_signals[1]))
                        {
                                if(ntk.is_on_critical_path(c) && !ntk.is_on_critical_path(d))
                                {
                                        //d = x3
                                        std::vector<signal> c_fanin_node_signals;
                                        node x, y;
                                        ntk.foreach_fanin(c, [&](auto const f)
                                        {
                                                c_fanin_node_signals.push_back(f);
                                        });
                                        x = ntk.get_node(c_fanin_node_signals[0]);
                                        y = ntk.get_node(c_fanin_node_signals[1]);

                                        if(ntk.is_on_critical_path(x) && !ntk.is_on_critical_path(y))
                                        {
                                                //y = x2, x = g
                                                //if(ntk.level(x) > (ntk.level(a) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(c_fanin_node_signals[1], n_fanin_node_signals[0]); //x2 x4
                                                        auto aig2 = ntk.create_and(c_fanin_node_signals[0], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!b_fanin_node_signals[1], n_fanin_node_signals[0]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth     
                                        }	
					else if(ntk.is_on_critical_path(y) && !ntk.is_on_critical_path(x))
                                        {
                                                //x = x2, y = g
                                                //if(ntk.level(y) > (ntk.level(a) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(c_fanin_node_signals[0], n_fanin_node_signals[0]); //x2 x4
                                                        auto aig2 = ntk.create_and(c_fanin_node_signals[1], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!b_fanin_node_signals[1], n_fanin_node_signals[0]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
                                        }
                                }
                                else if(ntk.is_on_critical_path(d) && !ntk.is_on_critical_path(c))
                                {
                                        //c = x3
                                        std::vector<signal> d_fanin_node_signals;
                                        node x, y;
                                        ntk.foreach_fanin(d, [&](auto const f)
                                        {
                                                d_fanin_node_signals.push_back(f);
                                        });
                                        x = ntk.get_node(d_fanin_node_signals[0]);
                                        y = ntk.get_node(d_fanin_node_signals[1]);

                                        if(ntk.is_on_critical_path(x) && !ntk.is_on_critical_path(y))
                                        {
                                                //y = x2, x = g
                                                //if(ntk.level(x) > (ntk.level(a) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(d_fanin_node_signals[1], n_fanin_node_signals[0]); //x2 x4
                                                        auto aig2 = ntk.create_and(d_fanin_node_signals[0], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!b_fanin_node_signals[0], n_fanin_node_signals[0]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
                                        }
					else if(ntk.is_on_critical_path(y) && !ntk.is_on_critical_path(x))
                                        {
                                                //x = x2, y = g
                                                //if(ntk.level(y) > (ntk.level(a) + 2))
                                                {
                                                        auto aig1 = ntk.create_and(d_fanin_node_signals[0], n_fanin_node_signals[0]); //x2 x4
                                                        auto aig2 = ntk.create_and(d_fanin_node_signals[1], aig1); // g (x2 x4)
                                                        auto aig3 = ntk.create_and(!b_fanin_node_signals[0], n_fanin_node_signals[0]); //x3 x4
                                                        auto aig4 = ntk.create_and(!aig2, !aig3);

                                                        ntk.substitute_node(n, !aig4); //g (x2 x4) + x3 x4
                                                	return true;
						}
                                                //else no benefit in terms of depth
                                        }
                                }
                        }
                        //else no distribution rule
		}

		return false;
	}
  }

private:
  Ntk& ntk;
};

} // namespace detail

/* Entry point for users to call */
template<class Ntk>
void aig_algebraic_rewriting( Ntk& ntk )
{
  static_assert( std::is_same_v<typename Ntk::base_type, aig_network>, "Ntk is not an AIG" );

  depth_view dntk{ntk};
  detail::aig_algebraic_rewriting_impl p( dntk );
  p.run();
}

} /* namespace mockturtle */
