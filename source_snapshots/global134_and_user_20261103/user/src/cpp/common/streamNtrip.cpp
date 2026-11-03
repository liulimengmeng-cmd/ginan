#include "common/zhangRealtimeNavigation.hpp"
#include <cstdlib>
// #pragma GCC optimize ("O0")

#include "common/streamNtrip.hpp"
#include <mutex>
#include "common/trace.hpp"

using std::lock_guard;
using std::mutex;

void TcpSocket::getData()
{
    lock_guard<mutex> guard(receivedDataBufferMtx);

    const int reserve = 8192;

    receivedData.reserve(reserve);

    for (auto it = chunkList.begin(); it != chunkList.end();)
    {
        auto& chunk = *it;

        if (receivedData.size() + chunk.size() > reserve)
        {
            break;
        }

        // 		std::cout << "\nCHUNK";
        // 		printHex(std::cout, *(vector<unsigned char>*)&chunk);

        receivedData.insert(receivedData.end(), chunk.begin(), chunk.end());

        it = chunkList.erase(it);
    }
}

void TcpSocket::dataChunkDownloaded(vector<char>& dataChunk)
{
    lock_guard<mutex> guard(receivedDataBufferMtx);

    if (zhangRealtimeNavigation::enabled())
    {
        std::size_t buffered = dataChunk.size() + receivedData.size();
        for (const auto& chunk : chunkList) buffered += chunk.size();
        if (dataChunk.size() > 8192 || buffered > 1024 * 1024)
        {
            std::cerr << "REALTIME_NTRIP_BUFFER_OVERFLOW generation_failed=1" << std::endl;
            std::_Exit(75);
        }
    }
    chunkList.push_back(std::move(dataChunk));
}
