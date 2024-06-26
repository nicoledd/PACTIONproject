/* Copyright 2004. Vladimir Prus
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE.txt or copy at
 * https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#include "../lists.h"
#include "../mem.h"
#include "../native.h"
#include "../object.h"
#include "../jam_strings.h"
#include "../variable.h"


/* Use quite klugy approach: when we add order dependency from 'a' to 'b', just
 * append 'b' to of value of variable 'a'.
 */
LIST * add_pair( FRAME * frame, int32_t flags )
{
    LIST * arg = lol_get( frame->args, 0 );
    LISTITER iter = list_begin( arg );
    LISTITER const end = list_end( arg );
    var_set( frame->module, list_item( iter ), list_copy_range( arg, list_next(
        iter ), end ), VAR_APPEND );
    return L0;
}


/* Given a list and a value, returns position of that value in the list, or -1
 * if not found.
 */
int32_t list_index( LIST * list, OBJECT * value )
{
    int32_t result = 0;
    LISTITER iter = list_begin( list );
    LISTITER const end = list_end( list );
    for ( ; iter != end; iter = list_next( iter ), ++result )
        if ( object_equal( list_item( iter ), value ) )
            return result;
    return -1;
}

enum colors { white, gray, black };


/* Main routine for topological sort. Calls itself recursively on all adjacent
 * vertices which were not yet visited. After that, 'current_vertex' is added to
 * '*result_ptr'.
 */
void do_ts( int32_t * * graph, int32_t current_vertex, int32_t * colors, int32_t * * result_ptr
    )
{
    int32_t i;

    colors[ current_vertex ] = gray;
    for ( i = 0; graph[ current_vertex ][ i ] != -1; ++i )
    {
        int32_t adjacent_vertex = graph[ current_vertex ][ i ];
        if ( colors[ adjacent_vertex ] == white )
            do_ts( graph, adjacent_vertex, colors, result_ptr );
        /* The vertex is either black, in which case we do not have to do
         * anything, or gray, in which case we have a loop. If we have a loop,
         * it is not clear what useful diagnostic we can emit, so we emit
         * nothing.
         */
    }
    colors[ current_vertex ] = black;
    **result_ptr = current_vertex;
    ( *result_ptr )++;
}


static void topological_sort( int32_t * * graph, int32_t num_vertices, int32_t * result )
{
    int32_t i;
    int32_t * colors = ( int32_t * )BJAM_CALLOC( num_vertices, sizeof( int32_t ) );
    for ( i = 0; i < num_vertices; ++i )
        colors[ i ] = white;

    for ( i = num_vertices - 1; i >= 0; --i )
        if ( colors[ i ] == white )
            do_ts( graph, i, colors, &result );

    BJAM_FREE( colors );
}


LIST * order( FRAME * frame, int32_t flags )
{
    LIST * arg = lol_get( frame->args, 0 );
    LIST * result = L0;
    int32_t src;
    LISTITER iter = list_begin( arg );
    LISTITER const end = list_end( arg );

    /* We need to create a graph of order dependencies between the passed
     * objects. We assume there are no duplicates passed to 'add_pair'.
     */
    int32_t length = list_length( arg );
    int32_t * * graph = ( int32_t * * )BJAM_CALLOC( length, sizeof( int32_t * ) );
    int32_t * order = ( int32_t * )BJAM_MALLOC( ( length + 1 ) * sizeof( int32_t ) );

    for ( src = 0; iter != end; iter = list_next( iter ), ++src )
    {
        /* For all objects this one depends upon, add elements to 'graph'. */
        LIST * dependencies = var_get( frame->module, list_item( iter ) );
        int32_t index = 0;
        LISTITER dep_iter = list_begin( dependencies );
        LISTITER const dep_end = list_end( dependencies );

        graph[ src ] = ( int32_t * )BJAM_CALLOC( list_length( dependencies ) + 1,
            sizeof( int32_t ) );
        for ( ; dep_iter != dep_end; dep_iter = list_next( dep_iter ) )
        {
            int32_t const dst = list_index( arg, list_item( dep_iter ) );
            if ( dst != -1 )
                graph[ src ][ index++ ] = dst;
        }
        graph[ src ][ index ] = -1;
    }

    topological_sort( graph, length, order );

    {
        int32_t index = length - 1;
        for ( ; index >= 0; --index )
        {
            int32_t i;
            LISTITER iter = list_begin( arg );
            for ( i = 0; i < order[ index ]; ++i, iter = list_next( iter ) );
            result = list_push_back( result, object_copy( list_item( iter ) ) );
        }
    }

    /* Clean up */
    {
        int32_t i;
        for ( i = 0; i < length; ++i )
            BJAM_FREE( graph[ i ] );
        BJAM_FREE( graph );
        BJAM_FREE( order );
    }

    return result;
}


void init_order()
{
    {
        char const * args[] = { "first", "second", 0 };
        declare_native_rule( "class@order", "add-pair", args, add_pair, 1 );
    }

    {
        char const * args[] = { "objects", "*", 0 };
        declare_native_rule( "class@order", "order", args, order, 1 );
    }
}
